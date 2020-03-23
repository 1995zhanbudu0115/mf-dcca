# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 15:13:31 2020

@author: zhanb
"""

import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings('ignore')

plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False


# 数据切片划分子序列
def devide(x,y,s,n,sort):
    xc_list = []
    yc_list = []
    if sort == 'desc':
        x = x[::-1].to_list()
        y = y[::-1].to_list()
    else:
        x = x.to_list()
        y = y.to_list()
    cut_index = [i for i in range(0, n, s)]
    pair_list = [(cut_index[i],cut_index[i+1]) for i in range(len(cut_index)-1)]
    for pair in pair_list:
        xc = x[pair[0]:pair[1]]
        yc = y[pair[0]:pair[1]]
        xc_list.append(xc)
        yc_list.append(yc)
    return xc_list,yc_list


# 一元n次拟合函数
def reg(x,y,order):
    coef = np.polyfit(x,y,order)
    poly_fit = np.polyval(coef,x)
    return poly_fit


# 对每个子序列拟合一元n次函数，并取拟合值
def all_ols(xv_list,i_list,order):
    xfit_list = []
    for xv in xv_list:
        xfit = list(reg(xv,i_list,order))
        xfit_list.append(xfit)
    return xfit_list


# 计算F2 前Ns个
def chige(arr):
    t_f2 = []
    for i in range(1, si+1):
        v1 = abs(arr[:si]*(arr[si]*si+i)-arr[si+1:si*2+1])
        v2 = abs(arr[si*2+1:si*3+1]*(arr[si]*si+i)-arr[si*3+1:])
        chi = np.dot(v1,v2)
        t_f2.append(chi)
    return np.average(t_f2)


# 计算F2 Ns+1到2Ns
def chige2(arr):
    t_f2 = []
    for i in range(1, si+1):
        v1 = abs(arr[:si]*(N-(arr[si]+1-Ns)*si+i)-arr[si+1:si*2+1])
        v2 = abs(arr[si*2+1:si*3+1]*(N-(arr[si]+1-Ns)*si+i)-arr[si*3+1:])
        chi = np.dot(v1,v2)
        t_f2.append(chi)
    return np.average(t_f2)


def linear_reg(x,y):
    linear = LinearRegression()
    linear.fit(x,y)
    coef = linear.coef_
    return coef


# 计算fq
def cfq(f2,q):
    if q != 0:
        fq = np.average((f2**(q/2)))**(1/q)
    else:
        fq = np.exp**(sum(f2)/2*len(f2))
    return fq


# 计算hq
def chq(fq,s):
    fq = np.array(fq)
    s = np.array(s)
    lns = np.log10(s).reshape(-1,1)
    lnfq = np.log10(fq).reshape(-1,1)
    hq = linear_reg(lns, lnfq)[0][0]
    print(hq)
    return hq


def all_fq(x_set,y_set,si):
    xi = x_set-x_set.mean()
    yi = y_set-y_set.mean()
    xv,yv = devide(xi,yi,si,N,sort='asd')
    xv2,yv2 = devide(xi,yi,si,N,sort='desc')
    i_list = [i for i in range(1,si+1)]
    xf_list = all_ols(xv,i_list,10)
    xf_list2 = all_ols(xv2,i_list,10)
    yf_list = all_ols(yv,i_list,10)
    yf_list2 = all_ols(yv2,i_list,10)
    xv_data = pd.DataFrame.from_records(xv)
    xv_data['v-1'] = xv_data.index
    xf_data = pd.DataFrame.from_records(xf_list)
    yv_data = pd.DataFrame.from_records(yv)
    yf_data = pd.DataFrame.from_records(yf_list)
    all_data = xv_data.merge(xf_data,left_index=True,right_index=True)
    all_data = all_data.merge(yv_data,left_index=True,right_index=True)
    all_data = all_data.merge(yf_data,left_index=True,right_index=True)
    all_values = all_data.values
    chi_arr = np.apply_along_axis(chige,1,all_values)
    xv_data2 = pd.DataFrame.from_records(xv2)
    xv_data2['v'] = xv_data.index+1+Ns
    xf_data2 = pd.DataFrame.from_records(xf_list2)
    yv_data2 = pd.DataFrame.from_records(yv2)
    yf_data2 = pd.DataFrame.from_records(yf_list2)
    all_data2 = xv_data2.merge(xf_data2,left_index=True,right_index=True)
    all_data2 = all_data2.merge(yv_data2,left_index=True,right_index=True)
    all_data2 = all_data2.merge(yf_data2,left_index=True,right_index=True)
    all_values2 = all_data2.values
    chi_arr2 = np.apply_along_axis(chige2,1,all_values2)
    all_chi = np.concatenate([chi_arr,chi_arr2])
    fq = cfq(all_chi,2)
    return fq


def rolling():
    data = pd.read_excel('data.xlsx')
    # 定义滚动参数
    t = 250
    values = data.iloc[:,[2,7]]
    hq_list = []
    # 进行滚动
    for i in range(len(data)-t):
        print(i,i+t)
        roll_data = values.iloc[i:i+t]
        x = roll_data.iloc[:,0]
        y = roll_data.iloc[:,1]
        global si
        fq_list = []
        si_list = [i for i in range(10, int(t/5))]
        # 循环 s ,s取值范围为10到t/5
        for si in range(10, int(t/5)):
            global N
            N = len(x)
            global Ns
            Ns = int(N/si)
            # 计算fq
            fq = all_fq(x,y,si)
            fq_list.append(fq)
        hq = chq(fq_list, si_list)
        hq_list.append(hq)
    hq_d = pd.DataFrame(hq_list, columns=['hq'], index=data.iloc[250:,0])
    hq_d.to_csv('hq.csv')
    hq_d.plot(figsize=(12,5))
    plt.savefig('hq.png')


if __name__ == '__main__':
    rolling()








    
    
                                                                                                                            