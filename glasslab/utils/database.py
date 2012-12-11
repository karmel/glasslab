'''
Created on Nov 15, 2010

@author: karmel
'''
from django.db import transaction, connections, connection
from subprocess import check_call
from glasslab.config import current_settings, django_settings
import datetime
from psycopg2 import OperationalError, Error as psycoError
import time

def execute_query(query, using='default', return_cursor=False):
    connection = connections[using]
    connection.close()
    cursor = connection.cursor()
    cursor.execute(query)
    transaction.commit_unless_managed()
    if return_cursor: return cursor
    connection.close()

def execute_query_without_transaction(query, using='default', return_cursor=False):
    connection = connections[using]
    isolation_level = connection.isolation_level
    connection.close()
    connection.isolation_level = 0
    cursor = connection.cursor()
    cursor.execute(query)
    transaction.commit_unless_managed()
    connection.isolation_level = isolation_level
    if return_cursor: return cursor
    connection.close()

def fetch_rows(query, return_cursor=False, using='default'):
    connection = connections[using]
    connection.close()
    cursor = connection.cursor()
    cursor.execute(query)
    result = cursor.fetchall()
    if return_cursor: return result, cursor
    connection.close()
    return result

def restart_server():
    '''
    On the Mac server, the Postgres server for some unknown reason
    does not release memory until restarted. We restart the server programmatically
    in certain circumstances to force the release of built of memory.
    
    To do this, we need to circumvent permissions issues by logging in 
    as user postgres. So, we ssh in with a public key set up on the 
    code-runner's machine (karmel@glass.bioinforma.tc, in this case)
    that is known to the postgres@glass.bioinforma.tc user. Then we restart 
    the server as the postgres user. 
    
    Notably, because we are restarting a daemon process, we have to 
    redirect all output into dev/null, or else the cursor doesn't return.
    
    We want to check first that the server knows it has restarted, and
    then that it can actually accept incoming queries.
    '''
    check_call('{0} "{1}/bin/pg_ctl restart -D {1}/data/ -m immediate </dev/null >/dev/null 2>&1 &"'.format(
                                            current_settings.PG_ACCESS_CMD, 
                                            current_settings.PG_HOME), shell=True)
    
    # Make sure everything went as planned
    time_to_wait = 5*60
    # We have to wait for the server to actually be up and running!
    server_is_starting = datetime.datetime.now()
    while server_is_starting:
        try:
            print fetch_rows('SELECT NOW();')
            connection.close()
            time.sleep(60)
            server_is_starting = False
        except OperationalError, psycoError:
            if (datetime.datetime.now() - server_is_starting).seconds > time_to_wait:
                server_is_starting = False
                raise Exception('Could not run query on server after {0} seconds. Please investigate.'.format(time_to_wait))
            
def discard_temp_tables(using='default'):
    '''
    Drops all temporary tables created in the current session.
    '''
    execute_query_without_transaction('DISCARD TEMP;', using=using)
    
    
    
class SqlGenerator(object):
    ''' 
    Parent class for schema-specific SQL generators.
    '''
    user = None
    def __init__(self, user=None):
        self.user = user or django_settings.DATABASES['default']['USER']
    
    def pkey_sequence_sql(self, schema_name, table_name):
        return """
        GRANT ALL ON TABLE "{0}"."{1}" TO  "{user}";
        CREATE SEQUENCE "{0}"."{1}_id_seq"
            START WITH 1
            INCREMENT BY 1
            NO MINVALUE
            NO MAXVALUE
            CACHE 1;
        ALTER SEQUENCE "{0}"."{1}_id_seq" OWNED BY "{0}"."{1}".id;
        ALTER TABLE "{0}"."{1}" ALTER COLUMN id SET DEFAULT nextval('"{0}"."{1}_id_seq"'::regclass);
        ALTER TABLE ONLY "{0}"."{1}" ADD CONSTRAINT {1}_pkey PRIMARY KEY (id);
        """.format(schema_name, table_name, user=self.user)