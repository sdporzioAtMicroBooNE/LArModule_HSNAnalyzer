#LArHSN - Liquid Argon Heavy Sterile Neutrinos

Development of LArSoft code for the analysis of heavy sterile neutrinos at the MicroBooNE detector.

In order to build it, you can use the same methods as for any other repository (e.g. LArSim).


After having a local install of LArSoft, you can just go to your `src` directory and then run:

```
mrb g https://github.com/sdporzioAtMicroBooNE/LArHSN
```

Then `cd` to your build directory (you can do `cd ${MRB_BUILDDIR}`).

```
mrbsetenv
mrb i
```

P.S.:
You can run `mrb uc` if you are having problems with your CMake files not updating.