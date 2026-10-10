import re
import settings

periodic_table = settings.periodic_table


def parse_formula(formula):
    elements = re.findall(r'([A-Z][a-z]*)(\d*\.\d+|\d*)', formula)
    element_dict = {}
    for element, amount in elements:
        amount = float(amount) if amount else 1.0
        element_dict[element] = amount
    return element_dict


def wt_to_at():
    while True:
        input_formula = input('请输入化学式（返回上一界面请按0）：')
        if input_formula == '0':
            break

        input_formula_dict = parse_formula(input_formula)

        print(f'你输入的元素成分为：{input_formula_dict}')

        n_list = [input_formula_dict[element] / periodic_table.get(element, 1) for element in input_formula_dict]
        n_all = sum(n_list)

        atom_and_at_percent_dict = {element: n / n_all * 100 for element, n in zip(input_formula_dict, n_list)}

        print(f'转化成at%为：{atom_and_at_percent_dict}')
        print(f'输入字典值的求和：{sum(input_formula_dict.values())}')
        print(f'输出字典值的求和：{sum(atom_and_at_percent_dict.values())}')


def at_to_wt():
    while True:
        input_formula = input('请输入化学式（返回上一界面请按0）：')
        if input_formula == '0':
            break

        input_formula_dict = parse_formula(input_formula)

        print(f'你输入的元素成分为：{input_formula_dict}')

        m_list = [input_formula_dict[element] * periodic_table.get(element, 1) for element in input_formula_dict]
        n_all = sum(m_list)

        atom_and_wt_percent_dict = {element: m / n_all * 100 for element, m in zip(input_formula_dict, m_list)}

        print(f'转化成wt%为：{atom_and_wt_percent_dict}')
        print(f'输入字典值的求和：{sum(input_formula_dict.values())}')
        print(f'输出字典值的求和：{sum(atom_and_wt_percent_dict.values())}')


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
