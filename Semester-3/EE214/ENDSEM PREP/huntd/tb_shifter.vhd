--------------------------------------------------------------
-- Test Bench for 3-bit shift register (ESD figure 2.6)
-- by Weijun Zhang, 04/2001
--
-- please note usually the processes within testbench do
-- not have sesitive list.
--------------------------------------------------------------	

library ieee;
use ieee.std_logic_1164.all;  
use ieee.std_logic_arith.all;
use ieee.std_logic_unsigned.all;

entity shifter_TB is			-- entity declaration
end shifter_TB;

architecture TB of shifter_TB is

    component shift_reg
    port(   I:	in std_logic;
	    clock:	in std_logic;
	    shift:	in std_logic;
	    Q:	out std_logic
    );
    end component;

    signal T_I:		std_logic;
    signal T_clock:	std_logic;
    signal T_shift:	std_logic;
    signal T_Q:		std_logic;

begin

    U_shifter: shift_reg port map (T_I, T_clock, T_shift, T_Q);
	
    -- concurrent process of clock
    process
    begin
	T_clock <= '0';
	wait for 5 ns;
	T_clock <= '1';
	wait for 5 ns;
    end process;

    -- concurrent process of test	
    process

	variable err_cnt: integer := 0;

    begin								
	T_shift <= '1';			-- start shifting  
	T_I <= '0';
	wait for 20 ns;
	T_I <= '1';	 		-- 1st/2nd bit input
	wait for 20 ns;
	T_I <= '0';			-- 3rd bit input
	wait for 10 ns;
	T_I <= '1';			-- 4th bit input
	wait;
    end process;
	
    process

	variable err_cnt: integer :=0; 

    begin
	
	-- case 1
	wait for 30 ns;
	assert(T_Q='0') report "Test1 Failed !"
	severity error;			 
	if (T_Q/='0') then
	    err_cnt:=err_cnt+1;
	end if;
		
	-- case 2
	wait for 10 ns;
	assert(T_Q='0') report "Test2 Failed !"
	severity error;
	if (T_Q/='0') then
	    err_cnt:=err_cnt+1;
	end if;
		
	-- case 3
	wait for 10 ns;
	assert(T_Q='1') report "Test3 Failed !"
	severity error;
	if (T_Q/='1') then
	    err_cnt:=err_cnt+1;
	end if;
		
	-- case 4
	wait for 10 ns;
	assert(T_Q='1') report "Test4 Failed !"
	severity error;
	if (T_Q/='1') then
	    err_cnt:=err_cnt+1;
	end if;
		
	-- summary of all the tests
	if (err_cnt=0) then 
	    assert (false) 
	    report "Testbench of Shifter completed successfully!"
	    severity note;
	else
	    assert (true)
	    report "Something wrong, try again!"
	    severity error;
	end if;
			
	wait;

    end process;

end TB;

----------------------------------------------------------------
configuration CFG_TB of shifter_TB is
	for TB
	end for;
end CFG_TB;
-----------------------------------------------------------------
