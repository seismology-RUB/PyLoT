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
        'wf':{
            'rgba': (0, 0, 0, 255)},
        'background': {
            'rgba': (255, 255, 255, 255)},
        'stylesheet': {
            'filename': None}
    },
    'dark':{
        'wf':{
            'rgba': (255, 255, 255, 255)},
        'background':{
            'rgba': (50, 50, 60, 255)},
        'stylesheet': {
            'filename': 'styles/dark.qss'}
    }
}

