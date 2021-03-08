module parser
    implicit none

    contains
        ! Get the file name provided via the command line
        subroutine get_config_file(filename)
            integer :: i
            character(len=32) :: arg
            character(:), allocatable, intent(out) :: filename

            filename = ''
            i = 0
            do
                call get_command_argument(i, arg)
                if (len_trim(arg) == 0) then
                    exit
                endif

                if (arg == '--config') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    filename = trim(arg)
                endif
                i = i+1
            end do

            if (filename == '') then
                write(*,*) 'No configuration file provided. Run code with "./epic --config [config file]"'
                return
            endif
        end subroutine get_config_file

end module parser
