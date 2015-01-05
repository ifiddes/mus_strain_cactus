import sqlite3 as sql


class ExclusiveSqlConnection(object):
    """meant to be used with a with statement to ensure proper closure"""
    def __init__(self, path):
        self.con = sql.connect(path, isolation_level = "EXCLUSIVE")
        con.execute("BEGIN EXCLUSIVE")
        self.cur = con.cursor()

    def __exit__(self):
        con.commit()
        con.close()


def hasTable(cur, table):
    """checks to make sure this sql database has a specific table"""
    cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='table_name'")
    rows = cur.fetchall()
    if table in rows:
        return True
    else:
        return False


def updateTable(cur, table, set_col_val, where):
    """
    table = name of table to be updated
    set_col_val = list of (column, value) pairs to be set
    where = pair of (column, value) where changes should occur
    """
    where_col, where_val = where
    for col, val in set_col_val:
        cur.execute("UPDATE {} SET {}=? WHERE {}=?".format(table, col, where_col),
                (val, where_val))


def initializeTable(cur, table, columns):
    """
    Builds a empty <table> in the sqlite db opened by <cur> with <[columns]>
    columns should be a list of name, type pairs
    """
    n, t = columns[0]
    cur.execute("CREATE TABLE {}({} {})".format(table, n, t))
    for n, t in columns[1:]:
        cur.execute("ALTER TABLE {} ADD COLUMN '{}' '{}'".format(table, n, t))
