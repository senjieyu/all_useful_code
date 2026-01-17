npt = '''variable dt equal 0.001
variable StepNPT equal 100000
variable PreStep equal 10000
variable nevery equal 100
variable out_dump_file_0 string force.0.dump
variable swap_attempt equal 300   # 每多少步尝试交换
variable swap_N equal 500        # 隔多少步交换
variable swap_temp equal 2000
# variable out_dump_file string
# variable mlip_ini string

units metal
boundary p p p
atom_style atomic

box tilt large
read_data data.in


pair_style  mlip ${mlip_ini}
pair_coeff * *

neighbor 2.0 bin
neigh_modify delay 10 check yes

timestep    ${dt}
variable Tdamp equal "v_dt * 100"
variable Pdamp equal "v_dt * 1000"
velocity        all create ${T} ${random} mom yes rot yes dist gaussian
thermo 5000


fix 100 all npt temp ${T} ${T} ${Tdamp} aniso 0.0 0.0 ${Pdamp}
velocity all scale ${T}
velocity all zero linear

run             ${PreStep}

reset_timestep 0

compute pe all pe/atom

group type1_atoms type 6
variable type1_count equal count(type1_atoms)
group type2_atoms type 8
variable type2_count equal count(type2_atoms)



if "${type1_count} > 0" then &
  "fix swap_5_6 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 5 6" &
  "fix swap_5_7 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 5 7" &
  "fix swap_5_9 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 5 9" &
  "fix swap_5_10 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 5 10" &
  "fix swap_6_7 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 6 7" &
  "fix swap_6_9 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 6 9" &
  "fix swap_6_10 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 6 10" &
  "fix swap_7_9 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 7 9" &
  "fix swap_7_10 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 7 10" &
  "fix swap_9_10 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 9 10" &
elif "${type2_count} > 0" &
  "fix swap_5_8 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 5 8" &
  "fix swap_5_7 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 5 7" &
  "fix swap_5_9 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 5 9" &
  "fix swap_5_10 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 5 10" &
  "fix swap_8_7 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 8 7" &
  "fix swap_8_9 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 8 9" &
  "fix swap_8_10 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 8 10" &
  "fix swap_7_9 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 7 9" &
  "fix swap_7_10 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 7 10" &
  "fix swap_9_10 all atom/swap ${swap_N} ${swap_attempt} ${random} ${swap_temp} ke yes types 9 10" &
else &
  "print 'NZSP'"


dump 0 all custom ${nevery} ${out_dump_file_0} id type x y z fx fy fz c_pe
dump_modify 0 format line "%6d %3d %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e" sort id

dump 1 all custom ${nevery} ${out_dump_file} id type x y z fx fy fz c_pe
dump_modify 1 sort id

run             ${StepNPT}
'''
