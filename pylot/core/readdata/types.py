#!/usr/bin/env python
#
# Provide user the opportunity to read arbitrary organized database
# types. This means e.g. seiscomp data structure (SDS) or event based
# EGELADOS structure. 
#

import os


class GenericDataBase(object):
	'''
	GenericDataBase type holds all information about the current data-
	base working on.
	'''
	def __init__(self, stexp=None, **kwargs):
		structExpression = os.path.split(stexp)
		self.dataBaseDict = kwargs 
		
				
		

class SeiscompDataStructure(GenericDataBase):
	pass
