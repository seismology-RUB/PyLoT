#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pylot.core.util.pdf import ProbabilityDensityFunction

pdf = ProbabilityDensityFunction.from_pick(0.34, 0.5, 0.54, type='exp')
pdf2 = ProbabilityDensityFunction.from_pick(0.34, 0.5, 0.54, type='exp')
diff = pdf - pdf2
