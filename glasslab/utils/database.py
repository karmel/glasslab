'''
Created on Oct 2, 2014

@author: karmel
'''
from pandas.io.sql import read_sql_query
from sqlalchemy import create_engine


def get_engine(uri, password, username='vespucci_user', dbname='vespucci'):
    '''
    dialect+driver://username:password@host:port/database
    '''
    return create_engine('postgresql://{}:{}@{}:5432/{}'.format(
        username, password, uri, dbname))


def dataframe_from_query(query, engine):
    '''
    Given the database connection and query, return a Pandas dataframe
    with the requested data.
    '''
    return read_sql_query(query, engine)
