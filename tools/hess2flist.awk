#!/usr/bin/gawk -f

#=========================================================================================
# Peacemaker -- A Quantum Cluster Equilibrium Code.
#
# Copyright 2004-2006 Barbara Kirchner, University of Bonn
# Copyright 2007-2012 Barbara Kirchner, University of Leipzig
# Copyright 2013-2016 Barbara Kirchner, University of Bonn
#
# This file is part of Peacemaker.
#
# Peacemaker is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Peacemaker is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Peacemaker.  If not, see <http://www.gnu.org/licenses/>
#=========================================================================================

/ir_spectrum/ {
    getline;
    n = $1;
    printf "%d\n\n", n;
    for (i = 0; i != n; ++i) {
        getline;
        printf "%f\n", $1;
    }
}
