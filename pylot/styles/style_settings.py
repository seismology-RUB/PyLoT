# -*- coding: utf-8 -*-

# Set base phase colors for manual and automatic picks
# together with a modifier (r, g, or b) used to alternate
# the base color
phasecolors = {
    'manual': {
        'P':{
            'rgba': (0, 0, 255, 255),
            'modifier': 'g'},
        'S':{
            'rgba': (255, 0, 0, 255),
            'modifier': 'b'}
    },
    'auto':{
        'P':{
            'rgba': (140, 0, 255, 255),
            'modifier': 'g'},
        'S':{
            'rgba': (255, 140, 0, 255),
            'modifier': 'b'}
    }
}

# Set plot colors and stylesheet for each style
stylecolors = {
    'default':{
        'linecolor':{
            'rgba': (0, 0, 0, 255)},
        'background': {
            'rgba': (255, 255, 255, 255)},
        'multicursor': {
            'rgba': (255, 190, 0, 255)},
        'ref': {
            'rgba': (200, 210, 230, 255)},
        'test': {
            'rgba': (200, 230, 200, 255)},
        'stylesheet': {
            'filename': None}
    },
    'dark': {
        'linecolor': {
            'rgba': (255, 255, 255, 255)},
        'background': {
            'rgba': (50, 50, 60, 255)},
        'multicursor': {
            'rgba': (0, 150, 190, 255)},
        'ref': {
            'rgba': (80, 110, 170, 255)},
        'test': {
            'rgba': (130, 190, 100, 255)},
        'stylesheet': {
            'filename': 'pylot/styles/dark.qss'}
    },
    'bright': {
        'linecolor': {
            'rgba': (0, 0, 0, 255)},
        'background': {
            'rgba': (255, 255, 255, 255)},
        'multicursor': {
            'rgba': (100, 100, 190, 255)},
        'ref': {
            'rgba': (200, 210, 230, 255)},
        'test': {
            'rgba': (200, 230, 200, 255)},
        'stylesheet': {
            'filename': 'pylot/styles/bright.qss'}
    }
}

