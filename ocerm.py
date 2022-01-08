# -*- coding: utf-8 -*-

'''
Author: Lorces
Pre-modify for OCERM files in Include
'''


import re
import os


ab_path = os.getcwd().replace('\\', '/')
AP = str(ab_path + '/CFD-Tools')


class IncludeCreate():
    # 创建初始入口边界条件文件
    def create_qbc(self, U, KB, H):
        qbc_o = open(str(AP + '/Include/Iteration_BC/infl.QBC'), 'w')
        visqbc_o = open(str(AP + '/Include/Iteration_BC/VIS.QBC'), 'w')

        line1 = '             1             1'
        line2 = '             2             1'

        q_per = str('%14.5f' % (100 / (KB - 1)))
        line1 = line1 + q_per * int(KB - 1)
        line2 = line2 + q_per * int(KB - 1)
        q_in = 0.95 * U * H
        q_in = str('{:>14.6f}'.format(q_in))
        vq_in = str('{:>14.6f}'.format(0))
        line3 = q_in * 2

        qbc_o.write(line1 + '\n')
        qbc_o.write(line2 + '\n')
        qbc_o.write('      0.000000\n')
        qbc_o.write(line3 + '\n')
        qbc_o.write('   9999.000000\n')
        qbc_o.write(line3 + '\n')
        print('infl.QBC for iteration has been createed.')

        ir = (KB - 1) % 8
        visqbc_o.write('      0.000000\n')
        for n in range(2):
            for i in range(1, KB, 8):
                if i <= (KB - 1 - ir):
                    line4 = 8 * vq_in
                else:
                    line4 = ir * vq_in
                visqbc_o.write(line4 + '\n')

        visqbc_o.write('   9999.000000\n')
        for n in range(2):
            for i in range(1, KB, 8):
                if i <= (KB - 1 - ir):
                    line4 = 8 * vq_in
                else:
                    line4 = ir * vq_in
                visqbc_o.write(line4 + '\n')

        print('VIS.QBC for iteration has been createed.')

        qbc_o.close()
        visqbc_o.close()

        print('Water depth H is: ' + str('{:.6f}'.format(H)) + ' and KB is ' + str(KB))
        return q_in

    # 生成监测点文件
    def create_gauge(self, NGAUGE, dx, x0, y0=0):
        gauge_o = open(str(AP + '/Output/Gauge_XY.dat'), 'w')
        gauge_o.write(str(NGAUGE) + '\n')

        for i in range(0, NGAUGE):
            x = x0 + dx * i
            y = y0
            gauge_o.write(str('{:8.4f}'.format(x)) + str('{:8.4f}'.format(y)) + '\n')

        print('Gauge_XY has been createed.')

    # 计算垂向分层
    def set_kb(self, KBH, H, Hdown, mult, mode=0):
        Hper = []
        dh = float(Hdown)
        Kmid, Hmid = 1, 1

        h = 0
        Hper.append(h)

        if mode == 0:
            for i in range(1, KBH):
                h += (H / (KBH - 1))
                Hper.append(h)
            KB = KBH
        else:
            for i in range(0, KBH):
                h += dh
                Hper.append(h)
                dh *= mult
            for i in range(1, 100):
                temp = (H - h) / i
                if temp < dh:
                    Kmid, Hmid = i - 1, h
                    break
            for i in range(0, Kmid):
                h += (H - Hmid) / Kmid
                Hper.append(h)
            KB = KBH + Kmid + 1

        for i in range(len(Hper)):
            Hper[i] = (1 - Hper[i] / H) * (-1)
        Hper.reverse()

        return KB, Hper


class IncludeModify():
    # 获取IJM,IJP,IJE
    def __init__(self):
        grd_i = open(str(AP + '/Include/OCERM.GRD'), 'r')
        cuv_i = open(str(AP + '/Include/OCERM.CUV'), 'r')

        ijm, ije, ijp = 0, 0, 0

        pat1 = 'ZONE'
        for line in grd_i:
            if re.search(pat1, line):
                ijp = int(re.search(r'I\s=\s*(\d*)', line).group(1))
                ijm = int(re.search(r'J\s=\s*(\d*)', line).group(1))

        pat2 = 'Detailed information of the edges'
        token = -1
        for index, line in enumerate(cuv_i):
            if re.search(pat2, line):
                token = index + 1
            if index == token:
                ije = int(re.search(r'\s*(\d*)', line).group(1))

        grd_i.close()
        cuv_i.close()

        self.IJM, self.IJE, self.IJP = ijm, ije, ijp

    # 修改网格几何信息文件OCERM.GRD
    def modify_grd(self, KB, H, Hper):
        grd_i = open(str(AP + '/Include/OCERM.GRD'), 'r')
        grd_o = open(str(AP + '/Output/OCERM.GRD'), 'w')

        grd_o.write('Grids and grid point depths and the topology of the meshes\n')
        grd_o.write('Vertical Segmentation_Sigma Levels (KB)\n')
        grd_o.write('   %d\n' % KB)

        for i in range(len(Hper)):
            grd_o.write('%8f\n' % Hper[i])

        grd_o.write('  Horizontal Segmentation  Levels IJP\n')
        grd_o.write('   ' + str(self.IJP) + '\n')

        pat1 = 'ZONE'
        token = 100
        for index, line in enumerate(grd_i):
            if re.search(pat1, line):
                token = index + 1
            if index >= token and index < (token + self.IJP):
                newline = line[0:32] + str('{:>16.6f}'.format(H)) + '\n'
                grd_o.write(newline)
            elif index >= (token + self.IJP):
                grd_o.write(line)

        grd_i.close()
        grd_o.close()

        print('OCERM.GRD has been modified.')

    # 修改网格几何信息文件OCERM.CUV
    def modify_cuv(self):
        cuv_i = open(str(AP + '/Include/OCERM.CUV'), 'r')
        cuv_o = open(str(AP + '/Output/OCERM.CUV'), 'w')

        for line in cuv_i:
            cuv_o.write(line)

        cuv_i.close()
        cuv_o.close()

        print('OCERM.CUV has been modified.')

    # 修改OCERM_INF头文件
    def modify_inf(self, KB, NGAUGE):
        inf_i = open(str(AP + '/Include/OCERM_INF'), 'r')
        inf_o = open(str(AP + '/Output/OCERM_INF'), 'w')

        pat1 = 'IJM='
        pat2 = 'NGAUGE='

        for line in inf_i:
            if re.search(pat1, line):
                line1 = '	  Parameter (IJM=%d,IJE=%d,IJP=%d,KB=%d,KSL=6)\n' % (self.IJM, self.IJE, self.IJP, KB)
                inf_o.write(line1)
            elif re.search(pat2, line):
                line2 = '	  Parameter (IVEL_INS=1,NGAUGE=%d)\n' % (NGAUGE)
                inf_o.write(line2)
            else:
                inf_o.write(line)

        inf_i.close()
        inf_o.close()

        print('OCERM_INF has been modified.')

    # 修改入口边界条件文件（结构化网格）
    def modify_qbc(self, q_in, KB):
        cuv_i = open(str(AP + '/Include/OCERM.CUV'), 'r')
        qbc_i = open(str(AP + '/Include/infl.QBC'), 'r')
        visqbc_i = open(str(AP + '/Include/VIS.QBC'), 'r')
        qbc_o = open(str(AP + '/Output/infl.QBC'), 'w')
        visqbc_o = open(str(AP + '/Output/VIS.QBC'), 'w')

        pat1 = '   Detailed information of the edges'
        pat2 = '-999'
        token = 10000000
        nx = 0

        for index, line in enumerate(cuv_i):
            if re.search(pat1, line):
                token = index + 1
            if token < index:
                if re.search(pat2, line) is None:
                    nx = index - token - 1
                    break

        line_qbc = qbc_i.readline()[14:]
        for i in range(1, nx + 1):
            line1 = '%14s' % str(i)
            qbc_o.write(line1 + line_qbc)

        qbc_o.write('      0.000000\n')
        ir = nx % 8
        for i in range(1, nx + 1, 8):
            line2 = '%14s' % str(q_in)
            if i <= (nx - ir):
                line2 = 8 * line2
            else:
                line2 = ir * line2
            qbc_o.write(line2 + '\n')

        qbc_o.write('   9999.000000\n')
        for i in range(1, nx + 1, 8):
            line2 = '%14s' % str(q_in)
            if i <= (nx - ir):
                line2 = 8 * line2
            else:
                line2 = ir * line2
            qbc_o.write(line2 + '\n')

        print('infl.QBC has been modified.')

        line_visqbc = ''
        if (KB - 1) % 8 == 0:
            upper = 2 * int((KB - 1) / 8)
        else:
            upper = 2 * (int((KB - 1) / 8) + 1)

        for index, line in enumerate(visqbc_i):
            if index >= 1 and index <= upper:
                line_visqbc = line_visqbc + line
            if index > upper:
                break

        visqbc_o.write('      0.000000\n')
        for i in range(1, nx + 1):
            visqbc_o.write(line_visqbc)

        visqbc_o.write('   9999.000000\n')
        for i in range(1, nx + 1):
            visqbc_o.write(line_visqbc)

        print('VIS.QBC has been modified.')

        cuv_i.close()
        qbc_i.close()
        visqbc_i.close()
        qbc_o.close()
        visqbc_o.close()

    # 修改出口边界条件文件（结构化网格）
    def modify_ebc(self, q_out):
        cuv_i = open(str(AP + '/Include/OCERM.CUV'), 'r')
        ebc_i = open(str(AP + '/Include/outl.EBC'), 'r')
        visebc_i = open(str(AP + '/Include/VIS.EBC'), 'r')
        ebc_o = open(str(AP + '/Output/outl.EBC'), 'w')
        visebc_o = open(str(AP + '/Output/VIS.EBC'), 'w')

        pat1 = '   Detailed information of the edges'
        pat2 = '\s+' + str(self.IJM) + '\s+' + str(self.IJM)
        pat3 = '\s+\d+\s+-999\s+\d+\s+\d+\s+\d+\s+-1'
        token1, token2 = 10000000, 10000000
        out_num = []

        for index, line in enumerate(cuv_i):
            if re.search(pat1, line):
                token1 = index + 1
            if token1 < index:
                if re.search(pat2, line):
                    token2 = index
            if token2 < index:
                if re.search(pat3, line) is None:
                    break
                num = int(re.search(r'\s+\d+\s+-999\s+(\d*)', line).group(1))
                out_num.append(num)

        for i in range(0, len(out_num)):
            line1 = '%14s' % str(out_num[i])
            ebc_o.write(line1 + '             3')
            if (i + 1) % 4 == 0 or i == (len(out_num) - 1):
                ebc_o.write('\n')

        nx = len(out_num)
        ebc_o.write('      0.000000\n')
        ir = nx % 8
        for i in range(1, nx + 1, 8):
            line2 = '%14s' % str(q_out)
            if i <= (nx - ir):
                line2 = 8 * line2
            else:
                line2 = ir * line2
            ebc_o.write(line2 + '\n')

        ebc_o.write('   9999.000000\n')
        for i in range(1, nx + 1, 8):
            line2 = '%14s' % str(q_out)
            if i <= (nx - ir):
                line2 = 8 * line2
            else:
                line2 = ir * line2
            ebc_o.write(line2 + '\n')

        print('outl.EBC has been modified.')

        line_visebc = '      0.000000'

        visebc_o.write('      0.000000\n')
        for j in range(2):
            for i in range(1, nx + 1):
                visebc_o.write(line_visebc)
                if i % 8 == 0 or i == nx:
                    visebc_o.write('\n')

        visebc_o.write('   9999.000000\n')
        for j in range(2):
            for i in range(1, nx + 1):
                visebc_o.write(line_visebc)
                if i % 8 == 0 or i == nx:
                    visebc_o.write('\n')

        print('VIS.EBC has been modified.')

        cuv_i.close()
        ebc_i.close()
        visebc_i.close()
        ebc_o.close()
        visebc_o.close()

    # 修改入口边界条件文件（非结构网格）
    def modify_qbc_unst(self, n_start, n_end, q_in, KB):
        cuv_i = open(str(AP + '/Include/OCERM.CUV'), 'r')
        qbc_i = open(str(AP + '/Include/infl.QBC'), 'r')
        visqbc_i = open(str(AP + '/Include/VIS.QBC'), 'r')
        qbc_o = open(str(AP + '/Output/infl.QBC'), 'w')
        visqbc_o = open(str(AP + '/Output/VIS.QBC'), 'w')

        line_qbc = qbc_i.readline()[14:]
        for i in range(n_start, n_end + 1):
            line1 = '%14s' % str(i)
            qbc_o.write(line1 + line_qbc)

        nx = n_end - n_start + 1

        qbc_o.write('      0.000000\n')
        ir = nx % 8
        for i in range(1, nx + 1, 8):
            line2 = '%14s' % str(q_in)
            if i <= (nx - ir):
                line2 = 8 * line2
            else:
                line2 = ir * line2
            qbc_o.write(line2 + '\n')

        qbc_o.write('   9999.000000\n')
        for i in range(1, nx + 1, 8):
            line2 = '%14s' % str(q_in)
            if i <= (nx - ir):
                line2 = 8 * line2
            else:
                line2 = ir * line2
            qbc_o.write(line2 + '\n')

        print('infl.QBC in unstructured mesh has been modified.')

        line_visqbc = ''
        if (KB - 1) % 8 == 0:
            upper = 2 * int((KB - 1) / 8)
        else:
            upper = 2 * (int((KB - 1) / 8) + 1)

        for index, line in enumerate(visqbc_i):
            if index >= 1 and index <= upper:
                line_visqbc = line_visqbc + line
            if index > upper:
                break

        visqbc_o.write('      0.000000\n')
        for i in range(1, nx + 1):
            visqbc_o.write(line_visqbc)

        visqbc_o.write('   9999.000000\n')
        for i in range(1, nx + 1):
            visqbc_o.write(line_visqbc)

        print('VIS.QBC in unstructured mesh has been modified.')

        cuv_i.close()
        qbc_i.close()
        visqbc_i.close()
        qbc_o.close()
        visqbc_o.close()

    # 修改出口边界条件文件（非结构网格）
    def modify_ebc_unst(self, n_start, n_end, q_out):
        cuv_i = open(str(AP + '/Include/OCERM.CUV'), 'r')
        ebc_i = open(str(AP + '/Include/outl.EBC'), 'r')
        visebc_i = open(str(AP + '/Include/VIS.EBC'), 'r')
        ebc_o = open(str(AP + '/Output/outl.EBC'), 'w')
        visebc_o = open(str(AP + '/Output/VIS.EBC'), 'w')

        for i in range(n_start, n_end + 1):
            line1 = '%14s' % str(i)
            ebc_o.write(line1 + '             1')
            if (i - n_start + 1) % 4 == 0 or i == (n_end):
                ebc_o.write('\n')

        nx = n_end - n_start + 1
        ebc_o.write('      0.000000\n')
        ir = nx % 8
        for i in range(1, nx + 1, 8):
            line2 = '%14s' % str(q_out)
            if i <= (nx - ir):
                line2 = 8 * line2
            else:
                line2 = ir * line2
            ebc_o.write(line2 + '\n')

        ebc_o.write('   9999.000000\n')
        for i in range(1, nx + 1, 8):
            line2 = '%14s' % str(q_out)
            if i <= (nx - ir):
                line2 = 8 * line2
            else:
                line2 = ir * line2
            ebc_o.write(line2 + '\n')

        print('outl.EBC in unstructured mesh has been modified.')

        line_visebc = '      0.000000'

        visebc_o.write('      0.000000\n')
        for j in range(2):
            for i in range(1, nx + 1):
                visebc_o.write(line_visebc)
                if i % 8 == 0 or i == nx:
                    visebc_o.write('\n')

        visebc_o.write('   9999.000000\n')
        for j in range(2):
            for i in range(1, nx + 1):
                visebc_o.write(line_visebc)
                if i % 8 == 0 or i == nx:
                    visebc_o.write('\n')

        print('VIS.EBC in unstructured mesh has been modified.')

        cuv_i.close()
        ebc_i.close()
        visebc_i.close()
        ebc_o.close()
        visebc_o.close()