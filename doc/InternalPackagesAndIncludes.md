
# Internal Packages and Includes

This memo specifies how internal dependencies are managed.
It does not address third party libraries (TPLs)

## General Rules
Do not use the special directory "..", and do not use angle brackets
expect for standard library includes, or perhaps system-wide includes.

## Packages
DMRG++ is composed of the following packages: PsimagLite, LanczosPlusPlus,
dmrg, and cincuenta. TODO: Maybe names should all be capitalized.

## How includes should be done
### Package includes one of its files
If package (say LanczosPlusPlus) wants to include a file
from within itself it should include it directly or may prepend it
with a directory internal to the package. These examples are correct.
```#include "MyFile.hpp"```
or
```#include "Engine/MyOtherFile.hpp"```

It must not use the special directory "..", and must not
prepend it with its own name.
These examples for package LanczosPlusPlus are all incorrect:
```#include "../MyFile.hpp"```
and
``#include "LanczosPlusPlus/MyFile.hpp"```
and
```#include "LanczosPlusPlus/Engine/MyOtherFile.hpp"```.

### Package includes a file from another Package
If a package wants to include a file from another package it must
prepend the include with package's name and inner location. But see the exception rule at the end.
These examples for package dmrg are correct:
``#include "LanczosPlusPlus/MyFile.hpp"```
or
``#include "LanczosPlusPlus/Engine/MyOtherFile.hpp"```

The examples for package dmrg are incorrect:
``#include "MyFile.hpp"```
and
``#include "Engine/MyOtherFile.hpp"```
and
``#include "../LanczosPlusPlus/MyFile.hpp"```

## Exceptions

### PsimagLite
The name PsimagLite can be dropped, and maybe Utilities can be dropped also.
These includes are allowed:
```#include "Utilities/Matrix.h"```
```#include "Matrix.h"```
```#include "IoNg.h"```


