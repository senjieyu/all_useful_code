import re
import settings

periodic_table = settings.periodic_table


def parse_formula(formula):
    """
    Parses a chemical formula and returns a dictionary with element symbols and corresponding percentages.
    Example: 'Bi27.51Cd60.07Pb12.42' => {'Bi': 27.51, 'Cd': 60.07, 'Pb': 12.42}
    """
    components = re.findall(r'([A-Z][a-z]*)(\d+\.\d+)', formula)
    return {element: float(percentage) for element, percentage in components}


def wt_to_at():
    print('你正在尝试将质量百分比转化为原子百分比')

    formula = input('请输入化学式：')
    input_atom_and_wt_percent_dic = parse_formula(formula)

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
    print(f'输入字典值的求和：{sum(input_atom_and_wt_percent_dic.values())}')
    print(f'输出字典值的求和：{sum(atom_and_at_percent_dic.values())}')


def at_to_wt():
    print('你正在尝试将原子百分比转化为质量百分比')

    formula = input('请输入化学式：')
    input_atom_and_at_percent_dic = parse_formula(formula)

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
    print(f'输入字典值的求和：{sum(input_atom_and_at_percent_dic.values())}')
    print(f'输出字典值的求和：{sum(atom_and_wt_percent_dic.values())}')


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
