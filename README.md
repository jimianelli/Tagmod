Tagmod
======

tm3.tpl is the basic model. Two batch files are included,
run.bat
runall.bat

runall with MCMC is invoked as   
`runall mc`
    
runall with just mode estimates is invoked as   
`runall `


Estimates are stored in archive directory

Revived 2018 model
------------------

The portable, validated one-area model is available in [`revived/`](revived/).
It builds with current ADMB releases on macOS and includes isolated run scripts,
a smoke test, the surviving 2018 inputs, archived parameter evidence, and
validation notes.

```sh
cd revived
./scripts/smoke-test.sh
```

The historical `tm3.tpl` and batch workflow remain at the repository root for
provenance.
