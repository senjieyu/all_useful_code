import settings

periodic_table = settings.periodic_table


def wt_to_at():
    # 成分输入，这里可以写成一个函数的形式，不过就两个逻辑，复制一下也无伤大雅
    input_atom_and_wt_percent_dic = {}
    print('你正在尝试将质量百分比转化为原子百分比')
    while True:
        input_atom = input('请输入元素（输入完毕请按0）：')
        if input_atom == '0':
            break
        elif input_atom in periodic_table:
            input_atom_wt = float(input('请输入元素对应的wt%：'))
            input_atom_and_wt_percent_dic[input_atom] = input_atom_wt
        else:
            print('请输入正确的元素化学式：')
            continue

    print(f'你输入的元素成分为：{input_atom_and_wt_percent_dic}')

    atom_list = []
    atom_wt_list = []
    for k, v in input_atom_and_wt_percent_dic.items():
        atom_list.append(k)
        atom_wt_list.append(v)

    n_list = []
    count = 0
    for i in atom_wt_list:
        n = i / periodic_table.get(atom_list[count])
        count += 1
        n_list.append(n)

    n_all = sum(n_list)

    atom_and_at_percent_dic = {}
    count = 0
    for i in atom_list:
        atom_and_at_percent_dic[i] = n_list[count] / n_all * 100
        count += 1

    print(f'转化成at%为：{atom_and_at_percent_dic}')


def at_to_wt():
    print('你正在尝试将原子百分比转化为质量百分比')
    input_atom_and_at_percent_dic = {}
    while True:
        input_atom = input('请输入元素（输入完毕请按0）：')
        if input_atom == '0':
            break
        elif input_atom in periodic_table:
            input_atom_at = float(input('请输入元素对应的at%：'))
            input_atom_and_at_percent_dic[input_atom] = input_atom_at
        else:
            print('请输入正确的元素化学式：')
            continue

    print(f'你输入的元素成分为：{input_atom_and_at_percent_dic}')

    atom_list = []
    atom_at_list = []
    for k, v in input_atom_and_at_percent_dic.items():
        atom_list.append(k)
        atom_at_list.append(v)

    m_list = []
    count = 0
    for i in atom_at_list:
        n = i * periodic_table.get(atom_list[count])
        count += 1
        m_list.append(n)

    n_all = sum(m_list)

    atom_and_wt_percent_dic = {}
    count = 0
    for i in atom_list:
        atom_and_wt_percent_dic[i] = m_list[count] / n_all * 100
        count += 1

    print(f'转化成at%为：{atom_and_wt_percent_dic}')


while True:
    choice = input('''
                1.将质量百分比转换成原子百分比wt%→at%
                2.将原子百分比转换成质量百分比at%→wt%
                3.退出
                请选择要进行的任务：''')
    if choice == '1':
        wt_to_at()
    elif choice == '2':
        at_to_wt()
    elif choice == '3':
        break
    else:
        print('请输入正确的号码！')
