from pathlib import Path

def merge_xyz_simple(folder: str = ".", out_name: str = "all_out_pure.xyz"):
    folder_path = Path(folder).resolve()
    out_path = folder_path / out_name

    # 找到当前目录下所有 .xyz（大小写不敏感），排除输出文件本身
    xyz_files = sorted(
        [p for p in folder_path.iterdir()
         if p.is_file()
         and p.suffix.lower() == ".xyz"
         and p.name != out_name],
        key=lambda p: p.name.lower()
    )

    print("目标目录:", folder_path)
    print("输出文件:", out_path)
    print("将合并这些文件：")
    for p in xyz_files:
        print("  -", p.name)

    if not xyz_files:
        print("没找到 .xyz 文件")
        return

    # 纯拼接：原样追加，不改任何内容
    with out_path.open("wb") as out:
        for p in xyz_files:
            out.write(p.read_bytes())

    print(f"✅ 完成：合并 {len(xyz_files)} 个文件")

if __name__ == "__main__":
    # 默认合并脚本所在目录（避免运行目录不对）
    script_dir = Path(__file__).resolve().parent
    merge_xyz_simple(str(script_dir))
