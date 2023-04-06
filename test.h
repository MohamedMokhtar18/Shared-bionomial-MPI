struct descr_t
{
    // Root (for collectives with root)
    int root = -1;

    // Buf size for current message (in bytes)
    int bufsize;

    // RMA window's id
    // win_id_t wid = MPI_WIN_NO_ID;
};
