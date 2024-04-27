-----------------------------------------------------------------
-- test bench for FSM (ESD book figure 2.7)
-- by Weijun Zhang, 04/2001
-----------------------------------------------------------------

library	ieee;
use ieee.std_logic_arith.all;  
use ieee.std_logic_unsigned.all;
use ieee.std_logic_1164.all;

entity fsm_tb is			-- entity declaration
end fsm_tb;

-----------------------------------------------------------------

architecture TB of fsm_tb is

signal T_a: std_logic;
signal T_clock: std_logic; 
signal T_reset:	std_logic;  
signal T_x: std_logic;

component seq_design
port(	a:		in std_logic;
	clock:		in std_logic;
	reset:		in std_logic;
	x:		out std_logic
);
end component;

begin
	
    U_fsm: seq_design port map(T_a, T_clock, T_reset, T_x);
	
    process
    begin
	T_clock <= '1';			-- clock cycle 10 ns
	wait for 5 ns;
	T_clock <= '0';
	wait for 5 ns;
    end process;
	
    process			   
	
	variable err_cnt: integer :=0;
	
    begin
	
	-- case 1		
	T_reset <= '1';
	wait for 20 ns;		
	assert (T_x='0') report "Failed Case 1" severity error;
	if (T_x/='0') then
	    err_cnt:=err_cnt+1;
	end if;
		
	-- case 2
	T_reset <= '0';
	T_a <= '0';											 
	wait for 20 ns;
	assert (T_x='0') report "Failed Case 2" severity error;
	if (T_x/='0') then
	    err_cnt:=err_cnt+1;
	end if;	
		
	-- case 3
	wait for 30 ns;
	T_a <= '1';
	wait for 35 ns;
	assert (T_x='1') report "Failed Case 3" severity error;
	if (T_x/='1') then
	    err_cnt:=err_cnt+1;
	end if;			
		
	-- case 4
	wait for 70 ns;
	T_reset <= '1';
	wait for 10 ns;
	assert (T_x='0') report "Failed Case 4" severity error; 
	if (T_x/='0') then
	    err_cnt:=err_cnt+1;
	end if;	
		
	-- summary of all the tests	
	if (err_cnt=0) then
	    assert false 
	    report "Testbench of FSM completely successfully!"
	    severity note;					
	else			
	    assert true
	    report "Something wrong, Check again pls!" 	
	    severity error;
	end if;
				
	wait;

    end process;	 
	
end TB;

----------------------------------------------------------------------
configuration CFG_TB of fsm_TB is
	for TB
	end for;
end CFG_TB;
-----------------------------------------------------------------------
