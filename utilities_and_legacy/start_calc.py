import os
start = 1
end = 7

pwd = os.getcwd()
jobs_script = os.path.join(pwd,'jobs_script')
os.chdir(jobs_script)
for i in range(start,end+1):
    #print(f'python bsub_{i}.lsf')
    if not os.path.exists(f'bsub_{i}.lsf'):
        print('errors')
    os.system(f'bsub<bsub_{i}.lsf')