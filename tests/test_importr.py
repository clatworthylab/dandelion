#!/usr/bin/env python

def test_importrpy2():
    try:
        from rpy2.robjects.packages import importr
        from rpy2.rinterface import NULL
        from rpy2.robjects import pandas2ri
    except:
        raise (ImportError(
            "Unable to initialise R instance. Please run this separately through R with Shazam's tutorial."
        ))

    sh = importr('shazam')
    assert sh.__module__ == 'rpy2.robjects.packages'
