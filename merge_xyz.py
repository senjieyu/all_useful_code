from pathlib import Path

def merge_xyz(folder=".", out_name="all.xyz"):
    folder_path = Path(folder).resolve()
    out_path = folder_path / out_name

    # 只选 .xyz / .extxyz（大小写不敏感），并排除输出文件
    xyz_files = sorted(
        [
            p for p in folder_path.iterdir()
            if p.is_file()
            and p.name != out_name
            and p.suffix.lower() in {".xyz", ".extxyz"}
        ],
        key=lambda p: p.name.lower()
    )

    print("目标目录:", folder_path)
    print("输出文件:", out_path)
    print("将合并这些文件：")
    for p in xyz_files:
        print("  -", p.name)

    if not xyz_files:
        print("没找到 .xyz/.extxyz 文件")
        return

    with out_path.open("wb") as out:
        for p in xyz_files:
            data = p.read_bytes()
            out.write(data)
            # 保险：确保文件之间至少有一个换行
            if not data.endswith(b"\n"):
                out.write(b"\n")

    print(f"✅ 完成：合并 {len(xyz_files)} 个文件 -> {out_path}")

if __name__ == "__main__":
    script_dir = Path(__file__).resolve().parent
    merge_xyz(str(script_dir), out_name="all.xyz")
