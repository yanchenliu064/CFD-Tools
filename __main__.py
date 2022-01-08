# -*- coding: utf-8 -*-

'''
Author: Lorces
Main program of CFD-Tools
'''

from mesh import *
from ocerm import *
from paperfig import *

# 1-Create mesh files, 2-Modify include files, 3-Create paper figure, 4-IBM visualization
module_name = 1

if module_name == 1:
    bd1 = BdCreate()
    boundary = bd1.add_arc([0, 0], 0.05, [0, 360], 100)

    rect1 = MeshFile('CYTE')
    rect1.rect_mesh(-1, 2, 61, -0.05, 0.05, 3)

    unst1 = MeshFile('CYM1')
    unst1.unst_mesh(boundary)

    custom1 = MeshCustom('CYM2')
    custom1.custom_mesh()

    ib1 = IbFile()
    ib1.ib_2D(boundary)
    ib1.ib_3D(3)

elif module_name == 2:
    u = 0.0004 * 0.05
    h = 0.1
    kb_down = 0.005
    kb_mode,kbh = 0,11
    kb_mult = 1.2

    gauge_num = 301
    gauge_x, gauge_y = -1, 2
    gauge_dx = 0.01

    inc_create = IncludeCreate()
    kb,per = inc_create.set_kb(kbh,h,kb_down,kb_mult,kb_mode)
    qIn = inc_create.create_qbc(u, kb, h)
    inc_create.create_gauge(gauge_num,gauge_dx,gauge_x)

    iteration = 0
    meshtype = 1

    ocerm1 = IncludeModify()

    if iteration == 0:
        ocerm1.modify_grd(kb,h,per)
        ocerm1.modify_cuv()
        ocerm1.modify_inf(kb,gauge_num)

        if meshtype == 1:
            ocerm1.modify_qbc(qIn,kb)
            ocerm1.modify_ebc('0.000000')
        elif meshtype == 2:
            in_start, in_end = 145,168
            out_start, out_end = 61,84
            ocerm1.modify_qbc_unst(in_start,in_end,qIn,kb)
            ocerm1.modify_ebc_unst(out_start,out_end,'0.000000')

elif module_name == 3:
    u_inf = 0.039
    d = 0.1
    gauge_begin = 1
    gauge_end = 160
    t_begin = 1000

    case1 = VisResult()

    case1.u_Cp()
    case1.u_mid(u_inf,d,gauge_begin,gauge_end,t_begin)
    case1.u_boundary()

elif module_name == 4:
    circle = [0, 0, 0.05]  # xr,yr,r,r2
    wmr = 0.012

    ibm_case1 = VisIBM()

    ibm_case1.ibm_gc_2D()
    ibm_case1.ibm_image_2D()
    ibm_case1.ibm_wm_2D(wmr, circle)
    ibm_case1.ibm_df_2D(circle)
    ibm_case1.ibm_ani_2D()

    ibm_case1.ibm_mesh(4)
    ibm_case1.ibm_vec()




