'''
Created on Sep 24, 2010

@author: karmel
'''

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'glasslab',
        'USER': 'glass',
        'PASSWORD': 'monocyte',
        'HOST': 'localhost',
        'PORT': '54321',
    },
    'read_only': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'glasslab',
        'USER': 'glass_read_only',
        'PASSWORD': 'monocyte',
        'HOST': 'localhost',
        'PORT': '54321',
    }
}
