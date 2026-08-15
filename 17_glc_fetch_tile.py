#!/usr/bin/env python3
"""
17_glc_fetch_tile.py — 从 Zenodo 上的 GLC_FCS30D zip 中只取出需要的瓦片

为什么要自己解 zip:
  GLC_FCS30D 在 Zenodo 按 5 度经度带打成 zip(共 194 GB)，中国范围要下约 50 GB。
  但每个 zip 内含 86 个 5x5 度瓦片(每个瓦片一个多波段 tif, 波段=2000..2022 各年)，
  站点实际只落在 46 个瓦片上，合计约 9 GB。
  GDAL 的 /vsizip//vsicurl/ 打不开(ZIP64 且 URL 以 /content 结尾, GDAL 认不出),
  所以这里自行解析中央目录、用 HTTP range 只下载需要的成员字节再解压。

用法:
  python3 17_glc_fetch_tile.py list  <E110>                # 列出该经度带包内的瓦片
  python3 17_glc_fetch_tile.py fetch <E110> <瓦片名> <输出tif>
"""
import sys, struct, zlib, os, time, subprocess

# 用 curl 而非 urllib: macOS 自带 Python 缺 CA bundle，urllib 会
# CERTIFICATE_VERIFY_FAILED；curl 用系统证书链，本项目全程可靠。

BASE = ("https://zenodo.org/api/records/8239305/files/"
        "GLC_FCS30D_19852022maps_{band}.zip/content")


def band_url(band):                       # band 形如 E110 -> E110-E115
    e = int(band[1:])
    return BASE.format(band=f"E{e}-E{e+5}")


def http_range(url, start, end, tries=5, exact=False):
    """取回 [start, end] 字节。

    注意不要给 curl 加 --retry: 断流重试时 curl 会把已输出的部分再输出一遍,
    字节流被悄悄地重复拼接, 解压时表现为 zlib "invalid block type" /
    "invalid code lengths set"。重试一律由这里整段重来。
    exact=True 时校验长度, 短读也算失败——解压是有状态的, 少一个字节就全废。
    """
    want = end - start + 1
    for k in range(tries):
        p = subprocess.run(["curl", "-sfL", "--max-time", "600",
                            "-r", f"{start}-{end}", url],
                           capture_output=True)
        if p.returncode == 0 and p.stdout and not (exact and len(p.stdout) != want):
            return p.stdout
        why = (f"长度 {len(p.stdout)} != {want}" if p.returncode == 0
               else f"curl {p.returncode}")
        print(f"    range {start}-{end} 第{k+1}次失败({why}), 重试", flush=True)
        time.sleep(3 * (k + 1))
    raise RuntimeError(f"range 请求失败: {start}-{end}")


def http_size(url):
    # Zenodo 的 HEAD 不给 Content-Length，用 1 字节 range 从 Content-Range 里读总长
    p = subprocess.run(["curl", "-sfL", "-D", "-", "-o", "/dev/null",
                        "--max-time", "120", "-r", "0-0", url],
                       capture_output=True, text=True)
    for ln in p.stdout.splitlines():
        if ln.lower().startswith("content-range:"):
            return int(ln.split("/")[-1].strip())
    raise RuntimeError("拿不到文件大小")


def central_dir(url):
    """返回 {成员名: (本地头偏移, 压缩大小, 压缩方法)}"""
    size = http_size(url)
    tail = http_range(url, max(0, size - 3_000_000), size - 1)

    i = tail.rfind(b"PK\x05\x06")
    if i < 0:
        raise RuntimeError("未找到 EOCD")
    cd_size, cd_off = struct.unpack("<II", tail[i + 12:i + 20])

    # ZIP64: 32 位字段被写成 0xFFFFFFFF 时，真实值在 ZIP64 EOCD 里
    if cd_off == 0xFFFFFFFF or cd_size == 0xFFFFFFFF:
        j = tail.rfind(b"PK\x06\x06")
        if j < 0:
            raise RuntimeError("未找到 ZIP64 EOCD")
        cd_size, cd_off = struct.unpack("<QQ", tail[j + 40:j + 56])

    cd = (tail[cd_off - (size - len(tail)):][:cd_size]
          if cd_off >= size - len(tail) else http_range(url, cd_off, cd_off + cd_size - 1))

    out, p = {}, 0
    while True:
        p = cd.find(b"PK\x01\x02", p)
        if p < 0:
            break
        method = struct.unpack("<H", cd[p + 10:p + 12])[0]
        csz, usz = struct.unpack("<II", cd[p + 20:p + 28])
        nlen, elen, clen = struct.unpack("<HHH", cd[p + 28:p + 34])
        lho = struct.unpack("<I", cd[p + 42:p + 46])[0]
        name = cd[p + 46:p + 46 + nlen].decode("utf-8", "replace")

        # 任一字段为 0xFFFFFFFF 说明真实值在 ZIP64 扩展字段里，按序补齐
        if 0xFFFFFFFF in (csz, usz, lho):
            ex = cd[p + 46 + nlen:p + 46 + nlen + elen]
            q = 0
            while q + 4 <= len(ex):
                hid, hsz = struct.unpack("<HH", ex[q:q + 4])
                if hid == 0x0001:
                    v, k = ex[q + 4:q + 4 + hsz], 0
                    if usz == 0xFFFFFFFF:
                        usz = struct.unpack("<Q", v[k:k + 8])[0]; k += 8
                    if csz == 0xFFFFFFFF:
                        csz = struct.unpack("<Q", v[k:k + 8])[0]; k += 8
                    if lho == 0xFFFFFFFF:
                        lho = struct.unpack("<Q", v[k:k + 8])[0]
                    break
                q += 4 + hsz
        out[name] = (lho, csz, method)
        p += 1
    return out


def fetch(url, name, lho, csz, method, dest):
    """下载单个成员的字节区间并解压"""
    head = http_range(url, lho, lho + 29)
    if head[:4] != b"PK\x03\x04":
        raise RuntimeError(f"本地头签名不对 @{lho}")
    nlen, elen = struct.unpack("<HH", head[26:30])
    data_off = lho + 30 + nlen + elen

    tmp = dest + ".part"
    dec = zlib.decompressobj(-15) if method == 8 else None
    got, CH = 0, 1 << 23
    with open(tmp, "wb") as f:
        while got < csz:
            end = min(got + CH, csz) - 1
            chunk = http_range(url, data_off + got, data_off + end, exact=True)
            got += len(chunk)
            f.write(dec.decompress(chunk) if dec else chunk)
            print(f"    {got/csz*100:5.1f}%  {got/1e6:.0f}/{csz/1e6:.0f} MB", flush=True)
        if dec:
            f.write(dec.flush())
    os.replace(tmp, dest)
    return os.path.getsize(dest)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    cmd, band = sys.argv[1], sys.argv[2]
    url = band_url(band)
    cdir = central_dir(url)
    if cmd == "list":
        for k in sorted(cdir):
            print(f"{k}\t{cdir[k][1]}")
    elif cmd == "fetch":
        name, dest = sys.argv[3], sys.argv[4]
        if name not in cdir:
            sys.exit(f"包内无此瓦片: {name}")
        lho, csz, method = cdir[name]
        print(f"  取 {name}  压缩 {csz/1e6:.0f} MB  方法 {method}", flush=True)
        n = fetch(url, name, lho, csz, method, dest)
        print(f"  完成 {dest}  {n/1e6:.0f} MB", flush=True)
