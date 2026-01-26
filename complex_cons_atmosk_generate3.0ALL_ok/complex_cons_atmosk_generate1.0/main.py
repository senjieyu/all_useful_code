# -*- coding: utf-8 -*-
# main.py

import sys
from setting import Config
from core_generator import run_generation

def main():
    print("==========================================")
    print("   Complex Twin Polycrystal Generator")
    print("==========================================")
    
    # 打印配置信息摘要
    print(f"Element: {Config.ELEMENT}")
    print(f"Twin Mode: {Config.TWIN_MODE}")
    print(f"Box Size: {Config.POLY_BOX}")
    print(f"Grains: {Config.POLY_GRAINS}")
    print(f"Output Dir: {Config.OUT_DIR}")
    print("------------------------------------------")
    
    try:
        run_generation(Config)
    except Exception as e:
        print(f"\n[ERROR] Generation failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
