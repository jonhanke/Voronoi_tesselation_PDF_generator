#*****************************************************************************
#       Copyright (C) 2013 Jonathan Hanke and Dan Yasaki
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#
##########################################################
## 
##  Generate a pdf picture of the well-rounded retract 
##  in a certain region of the complex upper-half plane
##  using NumPy and the PyX graphics library.
##  
##  Authors: Jonathan Hanke, Dan Yasaki -- 2013-09-18
##  
##########################################################

from pyx import *
from math import *
from numpy import *


def dist(a, b):
    """
    The distance between two points a and b in CC
    """
    c = a-b
    return sqrt((c * c.conjugate()).real)



def find_the_geodesic_arc_defined_by_two_points(a, b):
    """
    Take three points defining an arc in CC and return the arc
        a = start point
        b = end point
    """
    ## SPECIAL CASE: Deal with geodesic arcs that are not curving
    if abs(a.real - b.real) < 0.00001:
        return  path.line(a.real, a.imag,  b.real, b.imag)


    ## OTHERWISE: Find the center point in the real axis (or at infty)
    c = (a * a.conjugate() - b * b.conjugate()) / (2 * (a.real - b.real))

    ## Find the radius and centered endpoints
    r = dist(a,c)
    a1 = a-c
    b1 = b-c


    ## Only deal with upper angles!
    if a1.imag > 0:
        a_theta = acos((a1.real)/(1.0 * r)) * 180 / pi
    else:
        raise NotImplementedError
    if b1.imag > 0:
        try:
            b_theta = acos((b1.real)/(1.0 * r)) * 180 / pi
        except:
            raise RuntimeError, "Bad b1.real = " + str(b1.real) + ", r = " + str(r)
        
    else:
        raise NotImplementedError

    ## Create the shorter arc
    if a_theta <= b_theta:
        new_arc = path.path(path.arc(c.real, c.imag, r, a_theta, b_theta))
    else:
        new_arc = path.path(path.arc(c.real, c.imag, r, b_theta, a_theta))        

    ## Return the arc
    return new_arc



def generate_PSL2_translate_matrices_in_range(x_min, x_max, denom_max):
    """
    Gives the matrices [a, b; c, d] in SL_2(Z) with:
        x_min <= a/c, b/d <= x_max
        1 <= c,d <= denom_max
    or allowing c = 0.
    """
    M_list = []

    ## Generate the c = 0 case:
    for b in range(int(floor(x_min)), int(ceil(x_max))+1):
        M = matrix([[1, b], [0, 1]])
        M_list.append(M)

    ## Generate the c, d != 0 matrices
    for c in range(1, denom_max + 1):
        for d in range(1, denom_max + 1):
            for a in range(int(floor(x_min * c)), int(ceil(x_max * c)) + 1):
                for b in range(int(floor(x_min * d)), int(ceil(x_max * d)) + 1):
                    if a*d - b*c == 1:
                        M = matrix([[a, b], [c,d]])
                        M_list.append(M)

    ## Return the list of matrices
    return M_list


def generate_arcs_from_matrix_list(M_list):
    """
    List the arcs obtained by applying the given list of matrices 
    to the standard arc:

        find_the_arc_defined_by_three_points(rho, rho+1, 0)

    """
    arc_list = []
    rho = complex(-0.5, sqrt(3) * 0.5) 
    for M in M_list:
        a_new = (M[0,0] * rho + M[0,1]) / (M[1,0] * rho + M[1,1])
        b_new = (M[0,0] * (rho+1) + M[0,1]) / (M[1,0] * (rho+1) + M[1,1])
        new_arc = find_the_geodesic_arc_defined_by_two_points(a_new, b_new)
        arc_list.append(new_arc)

    ## Return the list of arcs
    return arc_list



## Set the parameters for our well-rounded retract picture
x_min = -0.6
x_max = 0.6
denom_max = 30

## Prepare the arcs (for the standard retract)
c = canvas.canvas()
M_list = generate_PSL2_translate_matrices_in_range(x_min, x_max, denom_max)
print "len(M_list) = ", len(M_list)
arc_list = generate_arcs_from_matrix_list(M_list)
print "len(arc_list) = ", len(arc_list)

## Generate and save the picture
for arc in arc_list:
    c.stroke(arc, [style.linewidth(0.005)])
c.writePDFfile("retract.pdf")


