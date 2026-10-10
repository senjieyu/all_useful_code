import os
import subprocess
import shutil
import sys

class GBGenerator:
    def __init__(self, atomsk_path=None):
        self.atomsk = atomsk_path or self._find_atomsk()
        if not self.atomsk:
            raise RuntimeError("Atomsk not found in PATH. Please ensure it is installed.")

    def _find_atomsk(self):
        # Try common paths or PATH
        path = shutil.which("atomsk")
        if path: return path
        
        # Check specific path from previous context
        specific_path = r"D:\2025.916aftersoftware\atmosk\Atomsk\atomsk.EXE"
        if os.path.exists(specific_path):
            return specific_path
            
        return None

    def run_cmd(self, cmd, cwd=None):
        print(f"[RUN] {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Atomsk failed: {e.stderr}")
            raise

    def _generate_core(self, input_file, output_folder, angle1, angle2, axis, stack_axis, overlap_dist, suffix="GB"):
        """
        Core logic: Rotate -> Rotate -> Merge -> Clean
        """
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        base_name = os.path.splitext(os.path.basename(input_file))[0]
        
        # 1. Rotate Grain 1
        # Convert to xsf for intermediate
        g1_file = os.path.join(output_folder, f"grain1_a{angle1}.xsf")
        cmd1 = [self.atomsk, input_file, "-rotate", axis, str(angle1), g1_file, "-overwrite"]
        self.run_cmd(cmd1)

        # 2. Rotate Grain 2
        g2_file = os.path.join(output_folder, f"grain2_a{angle2}.xsf")
        cmd2 = [self.atomsk, input_file, "-rotate", axis, str(angle2), g2_file, "-overwrite"]
        self.run_cmd(cmd2)

        # 3. Merge
        merged_file = os.path.join(output_folder, "merged_raw.xsf")
        # atomsk --merge <axis> 2 <file1> <file2> <output>
        cmd_merge = [self.atomsk, "--merge", stack_axis, "2", g1_file, g2_file, merged_file, "-overwrite"]
        self.run_cmd(cmd_merge)

        # 4. Remove Overlaps
        final_file = os.path.join(output_folder, f"{base_name}_{suffix}.lmp")
        # Using lmp as final output format
        cmd_clean = [self.atomsk, merged_file, "-remove-doubles", str(overlap_dist), final_file, "-overwrite"]
        self.run_cmd(cmd_clean)

        print(f"[DONE] Generated: {final_file}")
        return final_file

    def generate_tilt_gb(self, input_file, angle1, angle2, stack_axis="Y", tilt_axis="Z", overlap_dist=1.5):
        """
        Tilt GB: Rotation axis is perpendicular to stacking axis.
        Example: Stack Y, Rotate Z.
        """
        print(f"\n=== Generating Tilt GB (Stack: {stack_axis}, Axis: {tilt_axis}) ===")
        folder_name = f"Tilt_Stack{stack_axis}_Rot{tilt_axis}_A{angle1}_{angle2}"
        output_folder = os.path.join(os.path.dirname(input_file), "..", "output_GB", folder_name)
        
        return self._generate_core(input_file, output_folder, angle1, angle2, tilt_axis, stack_axis, overlap_dist, suffix="Tilt")

    def generate_twist_gb(self, input_file, angle1, angle2, stack_axis="Y", overlap_dist=1.5):
        """
        Twist GB: Rotation axis IS the stacking axis.
        Example: Stack Y, Rotate Y.
        """
        print(f"\n=== Generating Twist GB (Stack: {stack_axis}, Axis: {stack_axis}) ===")
        folder_name = f"Twist_Stack{stack_axis}_A{angle1}_{angle2}"
        output_folder = os.path.join(os.path.dirname(input_file), "..", "output_GB", folder_name)
        
        return self._generate_core(input_file, output_folder, angle1, angle2, stack_axis, stack_axis, overlap_dist, suffix="Twist")

    def generate_general_gb(self, input_file, angle1, angle2, rotation_axis, stack_axis="Y", overlap_dist=1.5):
        """
        General GB: Arbitrary rotation axis.
        rotation_axis can be "X", "Y", "Z" or "[hkl]" or "ux uy uz"
        """
        print(f"\n=== Generating General GB (Stack: {stack_axis}, Axis: {rotation_axis}) ===")
        # Sanitize axis string for folder name
        axis_str = str(rotation_axis).replace(" ", "").replace("[", "").replace("]", "")
        folder_name = f"General_Stack{stack_axis}_Rot{axis_str}_A{angle1}_{angle2}"
        output_folder = os.path.join(os.path.dirname(input_file), "..", "output_GB", folder_name)
        
        return self._generate_core(input_file, output_folder, angle1, angle2, rotation_axis, stack_axis, overlap_dist, suffix="General")

def main():
    # 示例用法
    base_dir = os.path.dirname(os.path.abspath(__file__))
    # 假设输入文件在 stru 目录下
    input_vasp = os.path.join(base_dir, "stru", "POSCAR_Cu111_twin_singleplane_d14_d210_d35_ortho.vasp")
    
    if not os.path.exists(input_vasp):
        print(f"Error: Input file not found: {input_vasp}")
        return

    gb_gen = GBGenerator()

    # 1. Tilt GB 示例: 绕 Z 轴旋转，沿 Y 轴堆垛
    # Grain 1: +10度, Grain 2: -10度
    gb_gen.generate_tilt_gb(input_vasp, angle1=10, angle2=-10, stack_axis="Y", tilt_axis="Z", overlap_dist=1.5)

    # 2. Twist GB 示例: 绕 Y 轴旋转，沿 Y 轴堆垛
    # Grain 1: 0度, Grain 2: 30度
    gb_gen.generate_twist_gb(input_vasp, angle1=0, angle2=30, stack_axis="Y", overlap_dist=1.5)

    # 3. General GB 示例: 绕 [111] 旋转 (假设为晶体坐标系), 沿 Y 轴堆垛
    # 注意: Atomsk -rotate 接受 Miller indices [hkl] 或 向量 x y z
    gb_gen.generate_general_gb(input_vasp, angle1=15, angle2=-15, rotation_axis="[111]", stack_axis="Y", overlap_dist=1.5)

if __name__ == "__main__":
    main()
