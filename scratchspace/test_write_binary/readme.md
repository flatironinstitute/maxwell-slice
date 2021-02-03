To run this test

First delete dd.bin (if exists)

Then generate dd.bin

```
./compile && ./test_write
```

Test reading it using:

test_read.m (in matlab)

or 

test_read.py (in python)

Note: there is no header on the fortran output, because we used `access='direct', recl=8`
