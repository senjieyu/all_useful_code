"""Double-clickable launcher for the CSL grain-boundary builder.

The editable JSON configuration is kept beside this script (or beside the
packaged .exe) in the ``configs`` folder.  This file also supports a small
command-line validation mode for automated checks.
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
import traceback
from dataclasses import asdict
from pathlib import Path
from tkinter import BOTH, END, LEFT, RIGHT, W, X, filedialog, messagebox, ttk
import tkinter as tk


IS_FROZEN = bool(getattr(sys, "frozen", False))
APP_ROOT = Path(sys.executable).resolve().parent if IS_FROZEN else Path(__file__).resolve().parent
SRC_ROOT = APP_ROOT / "src"
if not SRC_ROOT.is_dir() and not IS_FROZEN:
    raise RuntimeError(f"Source directory not found: {SRC_ROOT}")
if SRC_ROOT.is_dir():
    sys.path.insert(0, str(SRC_ROOT))

from csl_gb_builder.builder import build_grain_boundary  # noqa: E402
from csl_gb_builder.config import ConfigurationError, load_config  # noqa: E402


DEFAULT_CONFIG = APP_ROOT / "configs" / "sigma17_110_223_cu.json"


def validate_config(path: Path) -> str:
    """Load a configuration and return a concise displayable summary."""
    config = load_config(path)
    return json.dumps(asdict(config), ensure_ascii=False, indent=2, default=str)


class BuilderApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("CSL 晶界一键建模")
        self.minsize(760, 500)
        self.path_var = tk.StringVar(value=str(DEFAULT_CONFIG))
        self.status_var = tk.StringVar(value="选择或修改 configs 中的 JSON 配置后，点击“生成晶界”。")
        self._build_widgets()

    def _build_widgets(self) -> None:
        frame = ttk.Frame(self, padding=14)
        frame.pack(fill=BOTH, expand=True)
        ttk.Label(frame, text="配置文件（JSON）:").pack(anchor=W)
        row = ttk.Frame(frame)
        row.pack(fill=X, pady=(4, 10))
        ttk.Entry(row, textvariable=self.path_var).pack(side=LEFT, fill=X, expand=True)
        ttk.Button(row, text="选择…", command=self.choose_config).pack(side=RIGHT, padx=(8, 0))

        buttons = ttk.Frame(frame)
        buttons.pack(fill=X, pady=(0, 10))
        ttk.Button(buttons, text="校验配置", command=self.validate).pack(side=LEFT)
        self.generate_button = ttk.Button(buttons, text="生成晶界", command=self.generate)
        self.generate_button.pack(side=LEFT, padx=8)
        ttk.Button(buttons, text="打开配置文件夹", command=self.open_config_folder).pack(side=LEFT)

        ttk.Label(frame, textvariable=self.status_var, wraplength=720).pack(anchor=W, pady=(0, 8))
        self.output = tk.Text(frame, height=20, wrap="word")
        self.output.pack(fill=BOTH, expand=True)
        self.output.insert(END, "入口已就绪。默认示例配置：\n" + str(DEFAULT_CONFIG) + "\n")

    def _config_path(self) -> Path:
        return Path(self.path_var.get()).expanduser().resolve()

    def choose_config(self) -> None:
        filename = filedialog.askopenfilename(
            title="选择晶界配置", initialdir=str(APP_ROOT / "configs"),
            filetypes=[("JSON 配置", "*.json"), ("所有文件", "*.*")],
        )
        if filename:
            self.path_var.set(filename)

    def open_config_folder(self) -> None:
        folder = APP_ROOT / "configs"
        folder.mkdir(exist_ok=True)
        import os
        os.startfile(str(folder))

    def _show(self, content: str) -> None:
        self.output.delete("1.0", END)
        self.output.insert(END, content)

    def validate(self) -> None:
        try:
            path = self._config_path()
            self._show(validate_config(path))
            self.status_var.set("配置校验通过。")
        except (OSError, ConfigurationError, ValueError) as exc:
            self.status_var.set("配置校验失败。")
            self._show(str(exc))
            messagebox.showerror("配置错误", str(exc))

    def generate(self) -> None:
        try:
            config = load_config(self._config_path())
        except (OSError, ConfigurationError, ValueError) as exc:
            self.status_var.set("无法生成：配置无效。")
            messagebox.showerror("配置错误", str(exc))
            return
        self.generate_button.configure(state="disabled")
        self.status_var.set("正在调用 aimsgb 生成晶界，请稍候…")
        self._show("正在生成…\n")
        threading.Thread(target=self._generate_worker, args=(config,), daemon=True).start()

    def _generate_worker(self, config: object) -> None:
        try:
            result = build_grain_boundary(config)  # type: ignore[arg-type]
            text = json.dumps(result, ensure_ascii=False, indent=2, default=str)
            self.after(0, lambda: self._generation_done(text, None))
        except Exception as exc:  # Show dependency and backend failures to the user.
            detail = "".join(traceback.format_exception_only(type(exc), exc)).strip()
            self.after(0, lambda: self._generation_done(detail, exc))

    def _generation_done(self, text: str, error: Exception | None) -> None:
        self.generate_button.configure(state="normal")
        self._show(text)
        if error is None:
            self.status_var.set("生成完成；输出路径见上方 metadata。")
            messagebox.showinfo("完成", "晶界文件和 metadata 文件已生成。")
        else:
            self.status_var.set("生成失败。请确认 aimsgb 与 pymatgen 已按 README 安装。")
            messagebox.showerror("生成失败", text)


def main() -> int:
    parser = argparse.ArgumentParser(description="CSL grain-boundary builder launcher")
    parser.add_argument("--validate", metavar="CONFIG", type=Path, help="validate a JSON configuration and exit")
    args = parser.parse_args()
    if args.validate:
        try:
            print(validate_config(args.validate))
        except (OSError, ConfigurationError, ValueError) as exc:
            print(f"Configuration error: {exc}", file=sys.stderr)
            return 2
        return 0
    BuilderApp().mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
