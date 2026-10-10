'''银行相关业务的接口'''
from db import db_handle
def withdraw_interface(username,money):
    user_dic=db_handle.select(username)
    tot_money=money*1.05
    user_dic['balance'] -= tot_money
    if user_dic['balance'] >= 0 :
        #记录流水
        手续费 = tot_money - float(money)
        flow='取款成功,账户余额为{:.2f},手续费为{:.2f}'.format(user_dic['balance'],手续费)
        user_dic['flow'].append(flow)
        db_handle.save(user_dic)
        return True,flow
    else:
        return False, '提现金额超过现有账户余额，请重新输入'


def bank_repay(username,money):
    if money.isdigit() and money > '0' :
        money=float(money)
        user_dic=db_handle.select(username)
        user_dic['balance'] += money
        flow = '存款金额为{:.2f},账户余额为{:.2f}'.format(money,user_dic['balance'])
        user_dic['flow'].append(flow)
        db_handle.save(user_dic)
        return True,flow
    else:
        return False,'请输入正确的金额'



def bank_interface_transfer(username,transfer_user_name,money):
    from conf import settings
    # transfer_user_path = os.path.join(settings.BASE_PATH,'{}.json'.format(transfer_user_name))
    transfer_dic=db_handle.select(transfer_user_name)
    if transfer_dic:
        user_dic = db_handle.select(username)
        if user_dic['balance'] >= float(money):
            money=float(money)
            user_dic['balance'] -= money
            flow = "转账成功，转出金额：{:.2f},剩余金额：{:.2f}".format(money,user_dic['balance'])
            user_dic['flow'].append(flow)
            db_handle.save(user_dic)
            transfer_dic['balance'] += money
            flow1 = "{}向您转账，转入金额：{:.2f},剩余金额：{:.2f}".format(username,money, transfer_dic['balance'])
            transfer_dic['flow'].append(flow1)
            db_handle.save(transfer_dic)
            return True,flow
        else: return False,'余额不足'
    else:
        return False,'转入用户不存在'




