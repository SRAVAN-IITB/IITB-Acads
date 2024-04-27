library ieee;
use ieee.std_logic_1164.all;

entity system is
	port ( Clk_50 : in  std_logic;
		   Reset  : in  std_logic;
		   SW	  : in  std_logic;
		   LED	  : out std_logic_vector(7 downto 0));
end system;

architecture Struct of system is
	
	-- Do the component declaration if you need to use any other entity here before begin of architecture
	component Clock_Divider is
		port ( Clk_50  : in std_logic;
			   Reset   : in std_logic;
			   Clk_Out : out std_logic );
	end component Clock_Divider;
	
 	signal Clk_Out1s : std_logic;	
	begin
		-- MOST IMPORTANT NOTE: Component Instantantiation can not be done inside PROCESS block
		-- MOST IMPORTANT NOTE: Component Instantantiation can not be done inside PROCESS block
		-- MOST IMPORTANT NOTE: Component Instantantiation can not be done inside PROCESS block
		CLKDIV1 : Clock_Divider port map (Clk_50 => Clk_50, Reset => Reset, Clk_Out => Clk_Out); -- You can also control the reset signal for your desired
																								 -- functionality as you know when our system need to 
																								 -- restart counting
		
		L1 : process(Clk_Out)
				-- variable declaraction
				begin
					if (Clk_Out'event and Clk_Out = '1') then  -- Carefully notice you can the 250 Hz/500 Hz clock to generate delay of 1s / 0.5 sec
															   -- on the basis of theory how you created a 250Hz/500 Hz clock using 50 Mhz clock
						
						-- TASK 1 : Complete your desired functionality here to generate delay of 1s / 0.5s
						
						-- lets say Clk_Out1s is the signal generating delay of 1s/ 0.5s
					end if;
			 end process;
			 
		l2 : process (Clk_Out1s)
				-- Variable declaration
				begin
					if (Clk_Out1s'event and Clk_Out1s = '1') then
						if (SW = '0') then
							-- TASK 2 : Complete the part based on the fuctionality given to you for SW = 0
						else
							-- TASK 3 : Complete the part based on the fuctionality given to you for SW = 1
						end if;
 					end if;
					
					-- TASK 4 : Fill in the blanks based on the variable / signal you have used to generate the outputs of LED
					LED <= _____ ; 
			 end process;
			 
			 -- TASK 5 : Create your TESTBENCH to verify the following system
			 --			 To verify the system keep the delay less (not 1s / 0.5s) as it takes lot of time to observe the output in Modelsim 
			 --			 Use the exact given delay to test in your board
end Struct;