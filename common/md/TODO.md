# MDTrial Compatibility

* Add constructors
* Add convertors (reducers, really)
* Add method for mdplace!
* Make sure every generic method is truly generic
* Rename all references to MolecularDynamicsTrial to MDVerletTrial
* Abstract: 
    * `MDMetropolisTrial <: MDCanonicalTrial`
    * `MDVerletTrial <: MDMicroTrial <: MDCanonicalTrial
* Add abstract type tree:
  ```
  MDTrial
  └── MDCanonicalTrial
      ├── MDMetropolisTrial
      └── MDMicroTrial
          └── MDVerletTrial
  ```
