#coding: utf-8

import numpy as np

def simpleIntegLog(f, xlist):
    """
    当f为2维矩阵的时候矩阵的shape应为(n,m)
    其中n为xlist的长度
    """
    if hasattr(f, '__call__'):
        y = f(xlist)
    elif isinstance(f, np.ndarray):
        y = f
    else:
        raise RuntimeError('Error in simpleIntegLog!')
    y1 = y[:-1]
    y2 = y[1:]
    x1 = xlist[:-1]
    x2 = xlist[1:]
    if len(y1.shape) == 2:
        return (0.5*(y1.T*x1+y2.T*x2)*(np.log(x2)-np.log(x1))).sum(axis=1)
    return (0.5*(y1*x1+y2*x2)*(np.log(x2)-np.log(x1))).sum()

def folding(Fdnde,mymatrix,trueE,measureE='nouse'):
    """
    将真实能谱卷积响应矩阵
    Fdnde: 函数
    """
    if hasattr(Fdnde, '__call__'):
        y = Fdnde(trueE)
    elif isinstance(Fdnde, np.ndarray):
        y = Fdnde
    else:
        raise RuntimeError('Error in simpleIntegLog!')
    return simpleIntegLog((y*mymatrix.T).T,trueE)


def folding2(Fdnde,mymatrix,trueE,measureE='nouse'):
    """
    和folding是一样的
    用来作check
    """
    zz = []
    for i in range(len(trueE)):
        zz.append(Fdnde(trueE)[i]*mymatrix[i])
        
    zz = np.r_[zz].T
    dd1 = []
    for i,z_ in enumerate(zz):
        dd1.append(simpleIntegLog(z_,trueE))
    return np.r_[dd1]
