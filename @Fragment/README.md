Fragment
========

Template files
--------------

The template file specifies the geometry of the molecule. The atom labels have the elements as the first characters, the word ATOM and then the atom number afterwards.

You can also put "PAR#" with # starting at 1, anywhere in the file. This allows you to substitute values into this position in the file from inside a matlab script.

Example:

```text
 CATOM1
 OATOM2                  1      PAR1
 HATOM3                  1      1.11    2            122.2
 HATOM4                  1      1.11    2            122.2    3     - 180.0

!ENV
```

The file should be named *.tpl.

For a general basis, the format must be: (No newline between !ENV and the basis set, and the end must have 2 blank lines.

```text
CATOM1
 OATOM2                  1            1.23
 HATOM3                  1            1.11    2            122.2
 HATOM4                  1            1.11    2            122.2    3     -180.0

!ENV
-H     0
S   3   PAR1
      3.42525091             0.15432897
      0.62391373             0.53532814
      0.16885540             0.44463454
****
-C     0
S   3   1.00
     71.6168370              0.15432897
     13.0450960              0.53532814
      3.5305122              0.44463454
S   3   PAR2
      2.9412494             -0.09996723
      0.6834831              0.39951283
      0.2222899              0.70011547
P   3   PAR3
      2.9412494              0.15591627
      0.6834831              0.60768372
      0.2222899              0.39195739
****
-O     0
S   3   1.00
    130.7093200              0.15432897
     23.8088610              0.53532814
      6.4436083              0.44463454
S   3   PAR4
      5.0331513             -0.09996723
      1.1695961              0.39951283
      0.3803890              0.70011547
P   3   PAR5
      5.0331513              0.15591627
      1.1695961              0.60768372
      0.3803890              0.39195739
****


```

Create a fragment
-----------------

The matlab command:
    frag = Fragment('datapath', config);

creates a fragment, with the template and data stored in the 'datapath' directory. The config variable holds the specification for the calculation. The method `Fragment.defaultConfig()` returns a config structure with default values:

    template: 'template'  % Template file will be template.tpl.
    basisSet: 'STO-3G'    % Basis set keyword (Gaussian format).
    method: 'hf'          % Method to use ('hf' or 'mp2').
    charge: 0             % Charge on the fragment.
    spin: 1               % Spin (multiplicity) of the fragment, using Gaussian convention.
    par: []               % list of values to substitute for par1 par2 etc in the tpl file.

You then over-ride these anyway you want, i.e.

```matlab
config = Fragment.defaultConfig();
config.template = 'fhyde';
config.basisSet = 'sto-3g ';
config.par = [1.0 4.0];
frag = Fragment('c:\dave\apoly\msqc\data4', config);
```
