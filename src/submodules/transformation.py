
from abc import ABC


class TransformationManagerBase(ABC):
    '''
    TransformationManager abastract base class.

    Attributes
    ----------
    conn: sqlite.Connection
        A connection to a precursor database
    precursors: pd.DataFrame
        Long formatted dataframe of precursor quantities.
    proteins: pd.DataFrame
        Long formatted dataframe of protein quantities.
    '''

    def __init__(self, conn):
        self.conn = conn
        self.precursors = None
        self.proteins = None


    def get_long_tables(self, use_db_ids=False):
        '''
        Get long formatted precursor and protein tables.

        Parameters
        ----------
        use_db_ids: bool
            Should dataframe have database protein and replicate IDs?
        '''

        if use_db_ids:
            return self.precursors, self.proteins

        proteins = self.proteins.copy()
        precursors = self.precursors.copy()

        # add replicate column
        cur = self.conn.cursor()
        cur.execute('SELECT id, replicate FROM replicates;')
        rep_ids = {int(x[0]): x[1] for x in cur.fetchall()}
        precursors['replicate'] = precursors['replicateId'].apply(lambda x: rep_ids[x])
        proteins['replicate'] = proteins['replicateId'].apply(lambda x: rep_ids[x])

        # add protein name column
        cur.execute('SELECT proteinId, name FROM proteins;')
        prot_ids = {int(x[0]): x[1] for x in cur.fetchall()}
        proteins['protein'] = proteins['proteinId'].apply(lambda x: prot_ids[x])

        return precursors, proteins


    def get_wide_tables(self, normalized=True, use_db_ids=False):
        pass
