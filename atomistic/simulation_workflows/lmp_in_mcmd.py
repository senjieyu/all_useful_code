
# Original npt template (commented out)
# npt = '''variable dt equal 0.003
# variable StepNPT equal 100000
# variable PreStep equal 0
# variable nevery equal 100
# variable out_dump_file_0 string force.0.dump
# # variable out_dump_file string
# # variable mlip_ini string
#
# units metal
# boundary p p p
# atom_style atomic
#
# box tilt large
# read_data data.in
#
#
# pair_style  mlip ${mlip_ini}
# pair_coeff * *
#
# neighbor 2.0 bin
# neigh_modify delay 10 check yes
#
# timestep    ${dt}
# variable Tdamp equal "v_dt * 100"
# variable Pdamp equal "v_dt * 1000"
# velocity        all create ${T} ${random} mom yes rot yes dist gaussian
# thermo 100
#
# fix 100 all npt temp ${T} ${T} ${Tdamp} aniso 0.0 0.0 ${Pdamp}
# velocity all scale ${T}
# velocity all zero linear
#
# run             ${PreStep}
#
# reset_timestep 0
#
# compute pe all pe/atom
# dump 0 all custom ${nevery} ${out_dump_file_0} id type x y z fx fy fz c_pe
# dump_modify 0 format line "%6d %3d %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e" sort id
#
# dump 1 all custom ${nevery} ${out_dump_file} id type x y z fx fy fz c_pe
# dump_modify 1 sort id
#
# run             ${StepNPT}
# '''

# New npt template for Annealing (based on HfNbTiVAl_template.in)
npt = '''variable dt equal 0.001
variable StepNPT equal 8000000
variable StepEquil equal 40000
# T_start comes from external ${T}
variable T_start equal ${T}
variable T_end equal 300
variable P_target equal 0.0
variable nevery equal 4000
variable out_dump_file_0 string force.0.dump
# variable out_dump_file string
# variable mlip_ini string
# variable random equal 12345

units metal
boundary p p p
atom_style atomic


box tilt large
read_data data.in


# 【Core Fix 1】: Manually add masses for HfNbTiVAl (Solves ASE missing masses issue)
mass 1 26.9815 # Al
mass 2 47.867  # Ti
mass 3 50.9415 # V
mass 4 92.9064 # Nb
mass 5 178.49  # Hf



pair_style mlip ${mlip_ini}
pair_coeff * *

neighbor 2.0 bin
neigh_modify delay 10 check yes

timestep ${dt}
variable Tdamp equal "v_dt * 100"
variable Pdamp equal "v_dt * 1000"

# 【Core Fix 4】: Use rot no to avoid angular momentum issues with p p p boundary
velocity all create ${T_start} ${random} mom yes rot no dist gaussian

thermo ${nevery}
thermo_style custom step temp pe press vol lx ly lz

# ----------------- Equilibration Phase (Isothermal-Isobaric at T_start) -----------------
fix 100 all npt temp ${T_start} ${T_start} ${Tdamp} iso ${P_target} ${P_target} ${Pdamp}
run ${StepEquil}

reset_timestep 0

compute pe all pe/atom

# Dump 0 (forces reference)
dump 0 all custom ${nevery} ${out_dump_file_0} id type x y z fx fy fz c_pe
dump_modify 0 format line "%6d %3d %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e" sort id

# Dump 1 (trajectory)
dump 1 all custom ${nevery} ${out_dump_file} id type x y z fx fy fz c_pe
# 【Core Fix 5】: Element mapping for OVITO
dump_modify 1 element Al Ti V Nb Hf format float "%20.15e" sort id

# ----------------- Production Run: Annealing (T_start -> T_end) -----------------
unfix 100
fix 100 all npt temp ${T_start} ${T_end} ${Tdamp} iso ${P_target} ${P_target} ${Pdamp}
run ${StepNPT}

# Output final structure at 300K
write_data final_300K_structure.data
'''

# nvt = '''variable dt equal 0.003
# variable StepNVT equal 100000
# variable PreStep equal 0
# variable nevery equal 100
# variable out_dump_file_0 string force.0.dump
# # variable out_dump_file string
# # variable mlip_ini string
#
# units metal
# boundary p p p
# atom_style atomic
#
# box tilt large
# read_data data.in
#
#
# pair_style  mlip ${mlip_ini}
# pair_coeff * *
#
# neighbor 2.0 bin
# neigh_modify delay 10 check yes
#
# timestep    ${dt}
# variable Tdamp equal "v_dt * 100"
# variable Pdamp equal "v_dt * 1000"
# velocity        all create ${T} ${random} mom yes rot yes dist gaussian
# thermo 100
#
# fix 100 all nvt temp ${T} ${T} ${Tdamp}
# velocity all scale ${T}
# velocity all zero linear
#
# run             ${PreStep}
#
# reset_timestep 0
#
# compute pe all pe/atom
# dump 0 all custom ${nevery} ${out_dump_file_0} id type x y z fx fy fz c_pe
# dump_modify 0 format line "%6d %3d %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e" sort id
#
# dump 1 all custom ${nevery} ${out_dump_file} id type x y z fx fy fz c_pe
# dump_modify 1 sort id
#
# run             ${StepNVT}
# '''


nvt = '''variable dt equal 0.001
variable StepNVT equal 22000
variable PreStep equal 0
variable nevery equal 11
variable out_dump_file_0 string force.0.dump
# variable out_dump_file string
# variable mlip_ini string
# variable random equal 12345

units metal
boundary p p p
atom_style atomic

box tilt large
read_data data.in

# HfNbTiVAl masses
mass 1 26.9815
mass 2 47.867
mass 3 50.9415
mass 4 92.9064
mass 5 178.49

pair_style mlip ${mlip_ini}
pair_coeff * *

neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes

timestep ${dt}
variable Tdamp equal "v_dt * 100"

# 由外部接口传入温度和随机种子
velocity all create ${T} ${random} mom yes rot no dist gaussian

fix 100 all nvt temp ${T} ${T} ${Tdamp}

# 1步 swap + 10步 relax = 11步一个完整循环
thermo ${nevery}
thermo_style custom step temp pe etotal

# 如外部需要预平衡，可通过 PreStep 控制
run ${PreStep}

reset_timestep 0

compute pe all pe/atom

# 接口输出 0：参考 force dump
dump 0 all custom ${nevery} ${out_dump_file_0} id type x y z fx fy fz c_pe
dump_modify 0 format line "%6d %3d %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e" sort id

# 接口输出 1：主轨迹输出
dump 1 all custom ${nevery} ${out_dump_file} id type x y z fx fy fz c_pe
# 【Core Fix】: Element mapping for OVITO, matching npt
dump_modify 1 element Al Ti V Nb Hf format float "%20.15e" sort id

# ----------------- MC/MD hybrid loop -----------------
# 避免使用 loop 和 v_nframe 带来的浮点数解析错误
# 直接利用 fix atom/swap 的内置频率 N=11 实现每11步尝试1次swap
# 随机数种子使用 integer 格式确保 LAMMPS 不会报错

variable s12 equal ${random}+12
variable s12_int format s12 "%.0f"
variable s13 equal ${random}+13
variable s13_int format s13 "%.0f"
variable s14 equal ${random}+14
variable s14_int format s14 "%.0f"
variable s15 equal ${random}+15
variable s15_int format s15 "%.0f"
variable s23 equal ${random}+23
variable s23_int format s23 "%.0f"
variable s24 equal ${random}+24
variable s24_int format s24 "%.0f"
variable s25 equal ${random}+25
variable s25_int format s25 "%.0f"
variable s34 equal ${random}+34
variable s34_int format s34 "%.0f"
variable s35 equal ${random}+35
variable s35_int format s35 "%.0f"
variable s45 equal ${random}+45
variable s45_int format s45 "%.0f"

fix sw12 all atom/swap 11 1 ${s12_int} ${T} types 1 2
fix sw13 all atom/swap 11 1 ${s13_int} ${T} types 1 3
fix sw14 all atom/swap 11 1 ${s14_int} ${T} types 1 4
fix sw15 all atom/swap 11 1 ${s15_int} ${T} types 1 5
fix sw23 all atom/swap 11 1 ${s23_int} ${T} types 2 3
fix sw24 all atom/swap 11 1 ${s24_int} ${T} types 2 4
fix sw25 all atom/swap 11 1 ${s25_int} ${T} types 2 5
fix sw34 all atom/swap 11 1 ${s34_int} ${T} types 3 4
fix sw35 all atom/swap 11 1 ${s35_int} ${T} types 3 5
fix sw45 all atom/swap 11 1 ${s45_int} ${T} types 4 5

# 进行纯 NVT 并在指定步数自动进行 swap
run ${StepNVT}
'''
