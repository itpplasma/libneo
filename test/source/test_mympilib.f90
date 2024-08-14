program test_mympilib
  use derived_scheduler_module

  type(derived_scheduler) :: derived

  call derived%schedule()
end program test_mympilib
