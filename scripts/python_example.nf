process{
    input:
    file(my_csv) from csv_ch

    script:
    """
    #!/usr/bin/env python

    import pandas as pandas
    pd.read_csv($my_csv)
    """
}