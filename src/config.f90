!=========================================================================================
! Peacemaker -- A Quantum Cluster Equilibrium Code.
! 
! Copyright 2004-2006 Barbara Kirchner, University of Bonn
! Copyright 2007-2012 Barbara Kirchner, University of Leipzig
! Copyright 2013-2022 Barbara Kirchner, University of Bonn
!
! This file is part of Peacemaker.
! 
! Peacemaker is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! Peacemaker is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with Peacemaker.  If not, see <http://www.gnu.org/licenses/>
!=========================================================================================
! This module provides the config_t data type and associated functions, which are stored
! for parsing, storing, and accessing configuration files. A configuration file consists
! of sections, which hold records. Sections have a label, which must be unique. Records
! consist of keys and arguments. Keys must be unique within a section. Here is an example
! configuration file:
!
! # This is a comment
! [section1]
!   key1 arg1 arg2
!   key2
!   key3 arg1 arg2 arg3
!
! [section2]
!   # ...
!
! Records are stored in the record_t data type, whose public components provide access to
! the key, the number of arguments, and the arguments themselves.
!=========================================================================================
module config
    use iso_varying_string
    use error
    implicit none
    private
    !=====================================================================================
    ! Public entities.
    public :: config_t
    public :: record_t
    !=====================================================================================
    ! Represents a configuartion file.
    type :: config_t
        private
        type(section_t), pointer :: first_section => null()
        type(section_t), pointer :: last_section => null()
        integer :: n = 0
        contains
            procedure, public :: parse => config_read
            procedure, public :: has_section => config_has_section
            procedure, public :: check => config_check
            generic :: get_record => get_record_char, get_record_string
            procedure, private :: get_record_char => config_get_record_char
            procedure, private :: get_record_string => config_get_record_string
            procedure, public :: get_section_label => config_get_section_label
            procedure, public :: nr_sections => config_nr_sections
    end type config_t
    !=====================================================================================
    ! Represents a section.
    type :: section_t
        type(varying_string) :: label
        type(section_t), pointer :: next_section => null()
        type(record_t), pointer :: first_record => null()
        type(record_t), pointer :: last_record => null()
    end type section_t
    !=====================================================================================
    ! Represents a record.
    type :: record_t
        private
        integer, public :: nr_args
        type(varying_string), dimension(:), allocatable, public :: args
        type(varying_string), public :: key
        logical :: unread = .true.
        type(record_t), pointer :: next_record => null()
    end type record_t
    !=====================================================================================
    ! Character constants needed throughout this module.
    character, parameter :: tab = achar(9)
    character, parameter :: space = achar(32)
    character, parameter :: comment = '#' ! comment character
    character, parameter :: section_ld = '[' ! left section delimiter
    character, parameter :: section_rd = ']' ! right section delimiter
    !=====================================================================================
    contains
        !=================================================================================
        ! Reads a configuration file.
        subroutine config_read(config, filename)
            class(config_t), intent(inout) :: config
            character(*), intent(in) :: filename

            integer :: my_unit
            integer :: ios
            type(varying_string) :: line

            open(newunit = my_unit, file = filename, action = 'read', status = 'old', &
                iostat = ios)
            if (ios /= 0) call pmk_error("could not open '" // trim(filename) // "'")

            read_loop: do
                call get(my_unit, line, iostat = ios)
                if (ios == -1) exit read_loop
                call process_line(line)
                if (len(line) == 0) cycle read_loop
                if (is_section(line)) then
                    call new_section(config, line)
                else
                    call new_record(config, line)
                end if
            end do read_loop
        end subroutine config_read
        !=================================================================================
        ! Removes whitespace and comments. Converts tabs to spaces.
        subroutine process_line(line)
            type(varying_string), intent(inout) :: line

            integer :: comment_pos

            ! Replace tabs with spaces.
            line = replace(line, tab, space, every = .true.)
            ! Remove comments.
            comment_pos = index(line, comment)
            if (comment_pos /= 0) line = remove(line, comment_pos)
            ! Remove leading and trailing whitespaces.
            line = adjustl(line)
            line = trim(line)
        end subroutine process_line
        !=================================================================================
        ! Checks whether line is a section by scanning for left and right brackets.
        ! Performs syntax checks and stops if something is wrong.
        function is_section(line)
            type(varying_string), intent(in) :: line
            logical :: is_section

            integer :: ld_pos, rd_pos ! position of left and right section delimiter

            ! We scan for the position of the delimiters from the opposite end of the
            ! string. That way we can check for duplicate delimiters (syntax error) by
            ! checking whether the left delimiter is at the first position and
            ! whether the right delimiter is at the last position.
            ld_pos = index(line, section_ld, back = .true.)
            rd_pos = index(line, section_rd)

            is_section = .false.

            if (ld_pos == 0 .and. rd_pos == 0) return ! This is not a section.

            ! This looks like a section, check syntax.
            if (ld_pos == 1 .and. rd_pos == len(line)) then
                ! This is a section.
                is_section = .true.
            else
                ! There's a syntax error.
                call pmk_error("bad syntax in line '" // char(line) // "'")
            end if
        end function is_section
        !=================================================================================
        ! Adds a new section to the configuration. Section labels must be unique. Stops if
        ! a duplicate is encountered.
        subroutine new_section(config, line)
            type(config_t), intent(inout) :: config
            type(varying_string), intent(in) :: line

            type(section_t), pointer :: p
            type(varying_string) :: label

            ! Check for duplicate entries.
            label = extract(line, 2, len(line)-1)
            p => config%first_section
            do while (associated(p))
                if (p%label == label) call pmk_error("duplicate entry '[" // &
                    char(label) // "]'")
                p => p%next_section
            end do

            ! If we're here, there are no duplicates. Allocate space for a new section and 
            ! initialize it.
            if (.not. associated(config%first_section)) then
                allocate(config%first_section)
                config%last_section => config%first_section
            else
                allocate(config%last_section%next_section)
                config%last_section => config%last_section%next_section
            end if
            config%last_section%label = label
            config%n = config%n + 1
        end subroutine new_section
        !=================================================================================
        ! Checks whether there is a section with the given label within the configuration.
        function config_has_section(config, label)
            class(config_t), intent(in) :: config
            character(*), intent(in) :: label
            logical :: config_has_section

            type(section_t), pointer :: p

            p => config%first_section
            do while (associated(p))
                if (p%label == label) then
                    config_has_section = .true.
                    return
                end if
                p => p%next_section
            end do
            config_has_section = .false.
        end function config_has_section
        !=================================================================================
        ! Creates a new record.
        subroutine new_record(config, line)
            type(config_t), intent(inout) :: config
            type(varying_string), intent(inout) :: line

            type(section_t), pointer :: section ! pointer to current section
            type(record_t), pointer :: record ! pointer to current record
            type(varying_string) :: key
            integer :: i

            ! Check whether record is associated with any section.
            if (.not. associated(config%last_section)) call pmk_error("record '" // &
                char(line) // "is not associated with any section")

            ! Get the key and check that there is no duplicate record.
            call split(line, key, space)
            section => config%last_section
            record => section%first_record
            do while (associated(record))
                if (record%key == key) call pmk_error("duplicate entry '" // &
                    char(key) // "' in '[" // char(section%label) // "]'")
                record => record%next_record
            end do

            ! Allocate space for a new record_t data structure.
            if (.not. associated(section%first_record)) then
                allocate(section%first_record)
                section%last_record => section%first_record
            else
                allocate(section%last_record%next_record)
                section%last_record => section%last_record%next_record
            end if

            ! Split the remaining line and fill the record_t data structure.
            record => section%last_record
            record%key = key
            record%nr_args = count_words(line)
            allocate(record%args(record%nr_args))
            do i = 1, record%nr_args
                line = adjustl(line)
                call split(line, record%args(i), space)
            end do
        end subroutine new_record
        !=================================================================================
        ! Count the number of words in line.
        function count_words(line)
            type(varying_string), intent(in) :: line
            integer :: count_words

            character(len(line)) :: copy
            integer :: i
            logical :: in_word

            copy = line
            count_words = 0
            in_word = .false.
            do i = 1, len(copy)
                if (copy(i:i) == space) then
                    in_word = .false.
                else
                    if (.not. in_word) then
                        in_word = .true.
                        count_words = count_words + 1
                    end if
                end if
            end do
        end function count_words
        !=================================================================================
        ! Prints unread keys.
        subroutine config_check(config)
            class(config_t), intent(in) :: config

            type(section_t), pointer :: s
            type(record_t), pointer :: r

            s => config%first_section
            do while (associated(s))
                r => s%first_record
                do while (associated(r))
                    if (r%unread) call pmk_warning("unknown key '" // char(r%key) // &
                        "' in '[" // char(s%label) // "]'")
                    r => r%next_record
                end do
                s => s%next_section
            end do
        end subroutine config_check
        !=================================================================================
        ! Returns a pointer to the requested record, if it exists or null() if it does
        ! not. Sets the unread state of the record. A record is requested by a section
        ! label and a keyword.
        function config_get_record_char(config, section, key)
            class(config_t), intent(inout) :: config
            character(*), intent(in) :: section
            character(*), intent(in) :: key
            type(record_t), pointer :: config_get_record_char

            type(section_t), pointer :: s
            type(record_t), pointer :: r

            config_get_record_char => null()

            s => config%first_section
            do while (associated(s))
                if (s%label == section) then
                    r => s%first_record
                    do while (associated(r))
                        if (r%key == key) then
                            config_get_record_char => r
                            r%unread = .false.
                            return
                        end if
                        r => r%next_record
                    end do
                    ! Section labels are unique. It is pointless to continue searching.
                    return
                end if
                s => s%next_section
            end do
        end function config_get_record_char
        !=================================================================================
        ! Same as above, but accepts a string as section.
        function config_get_record_string(config, section, key)
            class(config_t), intent(inout) :: config
            type(varying_string), intent(in) :: section
            character(*), intent(in) :: key
            type(record_t), pointer :: config_get_record_string

            config_get_record_string => config_get_record_char(config, char(section), key)
        end function config_get_record_string
        !=================================================================================
        ! Returns the section label of the n'th section. Make sure that the section
        ! number is valid.
        function config_get_section_label(cfg, n)
            class(config_t), intent(in) :: cfg
            integer, intent(in) :: n
            type(varying_string) :: config_get_section_label

            integer :: i
            type(section_t), pointer :: s

            s => cfg%first_section
            do i = 2, n
                s => s%next_section
            end do
            config_get_section_label = s%label
        end function config_get_section_label
        !=================================================================================
        ! Returns the number of sections.
        function config_nr_sections(cfg)
            class(config_t), intent(in) :: cfg
            integer :: config_nr_sections

            config_nr_sections = cfg%n
        end function config_nr_sections
        !=================================================================================
end module config
!=========================================================================================
