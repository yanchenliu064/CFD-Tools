# -*- coding: utf-8 -*-

'''
Author: Lorces
Mesh file generation
'''


from math import *
import matplotlib.pyplot as plt
import numpy as np
import os

ab_path = os.getcwd().replace('\\', '/')
AP = str(ab_path + '/CFD-Tools')


class BdCreate():
    # 生成线段边界并剖分
    def add_line(self,l_start,l_end,n,is_show=None):
        x_start, y_start = l_start[0], l_start[1]
        x_end, y_end = l_end[0], l_end[1]
        dx = (x_end - x_start) / (n-1)
        dy = (y_end - y_start) / (n-1)
        points = np.zeros((n,2),dtype=float)

        for i in range(0,n):
            x = x_start + i*dx
            y = y_start + i*dy
            points[i,0] = x
            points[i,1] = y

        if is_show is not None:
            plt.figure()
            plt.plot(points[:,0],points[:,1],c='b',zorder=1)
            plt.scatter(points[:,0],points[:,1],c='r',zorder=2)
            plt.show()

        return points

    # 生成圆弧边界并剖分
    def add_arc(self,center,r,angle,n,is_show=None):
        xr, yr = center[0], center[1]
        angle_start, angle_end = angle[0], angle[1]
        angler = angle_end - angle_start
        points = np.zeros((n, 2), dtype=float)

        for i in range(0, n):
            a = (angle_start + (angler / n)) * i * 3.1415926 / 180
            x = xr + r * cos(a)
            y = yr + r * sin(a)
            points[i, 0] = x
            points[i, 1] = y

        if is_show is not None:
            plt.figure(figsize=(7,7))
            plt.plot(points[:,0],points[:,1],c='b',zorder=1)
            plt.scatter(points[:,0],points[:,1],c='r',zorder=2)
            plt.xlim((xr-1.1*r),(xr+1.1*r))
            plt.ylim((yr-1.1*r),(yr+1.1*r))
            plt.show()

        return points

    # 分段边界拼接
    def merge_boundary(self,boundary1,boundary2,is_end=None,is_show=None):
        if (boundary1[-1,0] == boundary2[0,0]) and (boundary1[-1,1] == boundary2[0,1]):
            boundary1 = np.delete(boundary1,-1,axis=0)

            if is_end is not None:
                boundary2 = np.delete(boundary2,-1,axis=0)

            boundary = np.concatenate([boundary1,boundary2])

            if is_show is not None:
                plt.figure(figsize=(7, 7))
                plt.plot(boundary[:, 0], boundary[:, 1], c='b', zorder=1)
                plt.scatter(boundary[:, 0], boundary[:, 1], c='r', zorder=2)
                plt.show()

            return boundary
        else:
            print('Boundaries cannot be merged, check it and try again.')

    # 边界简易可视化
    def final_boundary(self,boundary_out, boundary_in=None):
        plt.figure()
        plt.plot(boundary_out[:,0], boundary_out[:,1], c='b')
        plt.fill_between(boundary_out[:,0], boundary_out[:,1], facecolor='gray', alpha=0.3)
        plt.plot(boundary_in[:, 0], boundary_in[:, 1], c='b')
        plt.fill_between(boundary_in[:, 0], boundary_in[:, 1], facecolor='white')
        plt.show()


class MeshFile():
    def __init__(self,filename):
        self.filename = filename

    # 结构化网格文件生成
    def rect_mesh(self,xmin,xmax,nx,ymin,ymax,ny):
        rect = open(str(AP + '/Mesh/' + str(self.filename) + '.spl'), 'w')
        rect.write('%10s'% str(nx) + '%10s'%'2\n')
        xn = []
        h = (xmax - xmin) * 100 / (nx - 1)
        for i in range(0,nx):
            x = xmin + i*h/100
            xn.append(x)

        for i in range(len(xn)):
            line = '%10.4f' % xn[i] + '%10.4f' % ymin + '%10.4f' % xn[i] + '%10.4f' % ymax + '\n'
            rect.write(str(line))

        rect.write('%10s'% str(ny)+'\n')
        h = 1/(ny-1)
        for i in range(0,ny):
            y = i*h
            line = '%10.4f'% y + '\n'
            rect.write(line)

        print('Rect mesh '+str(self.filename)+' has been generated.')
        rect.close()

        xn, yn = [xmin,xmax,xmax,xmin,xmin], [ymin,ymin,ymax,ymax,ymin]
        plt.figure()
        plt.plot(xn,yn,c='b')
        plt.fill_between(xn, yn, facecolor='gray',alpha=0.3)
        plt.show()

    # 非结构网格文件生成
    def unst_mesh(self,boundary_out,boundary_in=None):
        unst = open(str(AP + '/Mesh/' + str(self.filename) +'.spl'), 'w')

        n_out = len(boundary_out)
        unst.write('%10s' % str(n_out) + '%10s' % '2\n')
        for i in range(0, n_out - 1):
            unst.write('%10.4f' % boundary_out[i,0] + '%10.4f' % boundary_out[i,1])
            unst.write('%10.4f' % boundary_out[i+1,0] + '%10.4f' % boundary_out[i+1,1] + '\n')
        unst.write('%10.4f' % boundary_out[n_out-1,0] + '%10.4f' % boundary_out[n_out-1,1])
        unst.write('%10.4f' % boundary_out[0,0] + '%10.4f' % boundary_out[0,1] + '\n')

        if boundary_in is not None:
            n_in = len(boundary_in)
            unst.write('%10s' % str(n_in) + '%10s' % '2\n')
            for i in range(0, n_in - 1):
                unst.write('%10.4f' % boundary_in[i, 0] + '%10.4f' % boundary_in[i, 1])
                unst.write('%10.4f' % boundary_in[i + 1, 0] + '%10.4f' % boundary_in[i + 1, 1] + '\n')
            unst.write('%10.4f' % boundary_in[n_in - 1, 0] + '%10.4f' % boundary_in[n_in - 1, 1])
            unst.write('%10.4f' % boundary_in[0, 0] + '%10.4f' % boundary_in[0, 1] + '\n')

        print('Unst mesh '+str(self.filename)+' has been generated.')

        unst.close()


class IbFile():
    # 二维虚拟离散边界数据文件生成
    def ib_2D(self,boundary):
        ib = open(str(AP + '/Mesh/IBM2D.dat'), 'w')

        n = len(boundary)

        ib.write('C1 \n')
        ib.write(str(n) + '\n')
        for i in range(0,n):
            ib.write('%10.5f' % boundary[i,0] + '%10.5f' % boundary[i,1] + '\n')

        plt.figure()
        plt.scatter(boundary[:,0], boundary[:,1], c='k')
        plt.show()

        print('2D immersed boundary has been generated.')
        ib.close()

    # 三维虚拟离散边界数据文件生成（ICFM-CFD文件格式修改）
    def ib_3D(self,npoint):
        index_raw = open(str(AP + '/Mesh/index_raw.dat'), 'r')
        index_new = open(str(AP + '/Mesh/index_new.dat'), 'w')
        ib = open(str(AP + '/Mesh/IBM3D.dat'), 'w')

        num = 0
        for line in index_raw:
            num = num + 1
            index_16 = line.strip('\n').split(' ')
            index = index_16.copy()
            for i in range(len(index_16)):
                index[i] = str(int(str(index_16[i]),16))
            index_w = ','.join(index)
            index_new.write(index_w+'\n')

        index_raw.close()
        index_new.close()

        num = num * 3
        index_r = open(str(AP + '/Mesh/index_new.dat'), 'r')
        pos = np.loadtxt(str(AP + '/Mesh/position.dat'))
        np.set_printoptions(formatter={'float': '{:16.8f}'.format})
        ib.write('%8s' % str(num) + '\n')

        for line in index_r:
            line = line.strip('\n')
            index = line.split(',')
            for i in range(npoint):
                pos_w = str(pos[int(index[i])-1,:])
                pos_w = ''.join(pos_w)
                ib.write(pos_w[1:-1] + '\n')

        print('3D immersed boundary has been modified from ICEM-CFD.')

        index_r.close()
        ib.close()


class MeshCustom():
    def __init__(self,filename):
        self.filename = filename

    # 自定义网格文件生成
    def custom_mesh(self):
        custom = open(str(AP + '/Mesh/' + str(self.filename) + '.spl'), 'w')

        xn1, yn1, xn2, yn2 = [], [], [], []

        for i in range(0, 60):
            x = -1 + i * 0.05
            y = -0.6
            xn1.append(x)
            yn1.append(y)

        for i in range(0, 24):
            x = 2
            y = -0.6 + i * 0.05
            xn1.append(x)
            yn1.append(y)

        for i in range(0, 60):
            x = 2 - i * 0.05
            y = 0.6
            xn1.append(x)
            yn1.append(y)

        for i in range(0, 24):
            x = -1
            y = 0.6 - i * 0.05
            xn1.append(x)
            yn1.append(y)

        custom.write('%10s' % str(len(xn1)) + '%10s' % '2\n')
        for i in range(0, len(xn1) - 1):
            custom.write('%10.4f' % xn1[i] + '%10.4f' % yn1[i])
            custom.write('%10.4f' % xn1[i + 1] + '%10.4f' % yn1[i + 1] + '\n')
        custom.write('%10.4f' % xn1[len(xn1) - 1] + '%10.4f' % yn1[len(xn1) - 1])
        custom.write('%10.4f' % xn1[0] + '%10.4f' % yn1[0] + '\n')


        xr, yr, r, nr = 0, 0, 0.16, 80
        for i in range(0, nr):
            a = - (360 / nr) * i * 3.1415926 / 180
            x = xr + r * cos(a)
            y = yr + r * sin(a)
            xn2.append(x)
            yn2.append(y)

        custom.write('%10s' % str(nr) + '%10s' % '2\n')
        for i in range(0, nr - 1):
            custom.write('%10.4f' % xn2[i] + '%10.4f' % yn2[i])
            custom.write('%10.4f' % xn2[i + 1] + '%10.4f' % yn2[i + 1] + '\n')
        custom.write('%10.4f' % xn2[nr - 1] + '%10.4f' % yn2[nr - 1])
        custom.write('%10.4f' % xn2[0] + '%10.4f' % yn2[0] + '\n')
        '''
    
        for i in range(0, 30):
            x = -0.15
            y = -0.15 + i * 0.01
            xn2.append(x)
            yn2.append(y)
    
        for i in range(0, 60):
            x = -0.15 + i * 0.01
            y = 0.15
            xn2.append(x)
            yn2.append(y)
    
        for i in range(0, 30):
            x = 0.45
            y = 0.15 - i * 0.01
            xn2.append(x)
            yn2.append(y)
    
        for i in range(0, 60):
            x = 0.45 - i * 0.01
            y = -0.15
            xn2.append(x)
            yn2.append(y)
    
        custom.write('%10s' % str(len(xn2)) + '%10s' % '2\n')
        for i in range(0, len(xn2) - 1):
            custom.write('%10.4f' % xn2[i] + '%10.4f' % yn2[i])
            custom.write('%10.4f' % xn2[i + 1] + '%10.4f' % yn2[i + 1] + '\n')
        custom.write('%10.4f' % xn2[len(xn2) - 1] + '%10.4f' % yn2[len(xn2) - 1])
        custom.write('%10.4f' % xn2[0] + '%10.4f' % yn2[0] + '\n')
        '''

        print('Custom mesh ' + str(self.filename) + ' has been generated.')
        custom.close()

        plt.figure()
        xn1.append(xn1[0])
        yn1.append(yn1[0])
        plt.plot(xn1, yn1, c='b')
        plt.fill_between(xn1, yn1, facecolor='gray', alpha=0.3)

        xn2.append(xn2[0])
        yn2.append(yn2[0])
        plt.plot(xn2, yn2, c='r')
        plt.fill_between(xn2, yn2, facecolor='white')
        plt.show()

    # 自定义虚拟边界文件生成
    def custom_ib(self):
        rec = np.zeros([2,5])
        rec[0, 0] = 0
        rec[1, 0] = -39
        rec[0, 1] = 344
        rec[1, 1] = -39
        rec[0, 2] = 344
        rec[1, 2] = 0
        rec[0, 3] = 0
        rec[1, 3] = 0
        rec[0, 4] = 9
        rec[1, 4] = 9

        theta = 156.2225 * np.pi / 180
        trsx = 5741.980273
        trsy = 12305.3496094

        rot = np.zeros([2,2])
        rot[0, 0] = np.cos(theta)
        rot[0, 1] = -1 * np.sin(theta)
        rot[1, 0] = np.sin(theta)
        rot[1, 1] = np.cos(theta)

        ibm = np.dot(rot,rec)

        ibm[0, :] = ibm[0, :] + trsx
        ibm[1, :] = ibm[1, :] + trsy

        ib = open(str(AP + '/Mesh/IBM2D.dat'), 'w')

        for i in range(5):
            ib.write('%16.6f' % ibm[0,i] + '\t'+'%16.6f' % ibm[1,i] + '\n')

    # 柱面网格文件生成
    def cycl_mesh(self):
        cycl = open(str(AP + '/Mesh/' + str(self.filename) + '.spl'), 'w')

        nx = 80

        cycl.write('%10s'% str(nx)+'\n')
        xn1, yn1 = [], []
        xn2, yn2 = [], []

        nr = 80
        xr = 0
        yr = 0
        r = 0.04
        for i in range(0, nr):
            a = -1 * (360 / nr) * i * 3.1415926 / 180
            x = xr + r * cos(a)
            y = yr + r * sin(a)
            xn1.append(x)
            yn1.append(y)

        nr = 80
        xr = 0
        yr = 0
        r = 0.16
        for i in range(0, nr):
            a = -1 * (360 / nr) * i * 3.1415926 / 180
            x = xr + r * cos(a)
            y = yr + r * sin(a)
            xn2.append(x)
            yn2.append(y)


        for i in range(0, len(xn1)):
            cycl.write('%10.4f' % xn1[i] + '%10.4f' % yn1[i])
            cycl.write('%10.4f' % xn2[i] + '%10.4f' % yn2[i] + '\n')

        ny = 21
        cycl.write('%10s'% str(ny)+'\n')
        h = 1/(ny-1)
        for i in range(0,ny):
            y = i*h
            line = '%10.4f'% y + '\n'
            cycl.write(line)

        print('Cycl mesh '+str(self.filename)+' has been generated.')
        cycl.close()

    # 局部加密网格文件生成
    def refine_mesh(self):
        refine = open(str(AP + '/Mesh/' + str(self.filename) + '.spl'), 'w')

        xn, xn1, xn2, xn3, xn4, xn5 = [], [], [], [], [], []
        yn, yn1, yn2, yn3, yn4, yn5 = [], [], [], [], [], []
        un, un1, un2, un3, un4, un5 = [], [], [], [], [], []

        h = 0.1
        uin = 0.00045
        ymin = -0.8
        ymax = 0.8

        dmin = 0.005
        dnum = 10
        dh = 2 * (0.2-dmin*dnum) / (dnum*(dnum+1))

        for i in range(0,int(0.6/dmin+1)):
            x = -0.3 + i * dmin
            y = -0.3 + i * dmin
            u = 0.95 * uin * h * dmin
            xn1.append(x)
            yn1.append(y)
            un1.append(u)
        un1.pop()

        x = -0.3
        y = -0.3
        for i in range(1,int(dnum+1)):
            dx = dmin + i * dh
            x = x - dx
            dy = dmin + i * dh
            y = y - dy
            u = 0.95 * uin * h * dy
            xn2.append(x)
            yn2.append(y)
            un2.append(u)
        xn2.reverse()
        yn2.reverse()
        un2.reverse()

        x = 0.3
        y = 0.3
        for i in range(1, int(dnum+1)):
            dx = dmin + i * dh
            x = x + dx
            dy = dmin + i * dh
            y = y + dy
            u = 0.95 * uin * h * dy
            xn3.append(x)
            yn3.append(y)
            un3.append(u)

        for i in range(1, 11):
            x = -0.5 - i*0.05
            xn4.append(x)
        xn4.reverse()

        for i in range(1, 31):
            x = 0.5 + i*0.05
            xn5.append(x)

        for i in range(1, 7):
            y = -0.5 - i * 0.05
            u = 0.95 * uin * h * 0.05
            yn4.append(y)
            un4.append(u)
        yn4.reverse()
        un4.reverse()

        for i in range(1, 7):
            y = 0.5 + i * 0.05
            u = 0.95 * uin * h * 0.05
            yn5.append(y)
            un5.append(u)

        xn = xn4 + xn2 + xn1 + xn3 + xn5
        yn = yn4 + yn2 + yn1 + yn3 + yn5
        un = un4 + un2 + un1 + un3 + un5

        nx = len(xn)
        ny = len(yn)

        refine.write('%10s' % str(nx) + '%10s' % '2\n')
        for i in range(len(xn)):
            line = '%10.4f' % xn[i] + '%10.4f' % ymin + '%10.4f' % xn[i] + '%10.4f' % ymax + '\n'
            refine.write(str(line))

        refine.write('%10s' % str(ny) + '\n')
        for i in range(0,len(yn)):
            posy = yn[i]/(ymax-ymin)+0.5
            line = '%10.4f'% posy + '\n'
            refine.write(line)
        udata = open(str(AP + '/Mesh/udata.dat'), 'w')

        for i in range(len(un)):
            line = '{:>16.8f}'.format(un[i])
            udata.write(line)
            if ((i+1) % 8 == 0):
                udata.write('\n')
        print('Rect mesh '+str(self.filename)+' has been generated.')

        refine.close()


