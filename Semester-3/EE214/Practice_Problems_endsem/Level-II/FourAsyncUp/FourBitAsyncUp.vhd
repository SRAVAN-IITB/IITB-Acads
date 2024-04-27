library IEEE;
use IEEE.std_logic_1164.all;

entity FourBitAsyncUp is
	port (Clock : in std_logic;
			Counter : out std_logic_vector(3 downto 0);
			Enable, Preset, Reset : in std_logic);
end entity FourBitAsyncUp;

architecture Struct of FourBitAsyncUp is
	component JK_FlipFlop is
	port (J, K : in std_logic;
			Clock : in std_logic;
			Enable : in std_logic;
			Preset, Reset : in std_logic;
			Q, NotQ : out std_logic);
	end component JK_FlipFlop;
	
	signal Q : std_logic_vector(3 downto 0) := "0000";
	signal NotQ : std_logic_vector(3 downto 0) := "1111";
	
	begin
	JKFF1 : JK_FlipFlop port map (J => '1', K => '1',
											Clock => Clock,
											Enable => Enable,
											Preset => Preset,
											Reset => Reset,
											Q => Q(0), NotQ => NotQ(0));
											
	GenFF : for i in 1 to 3 generate
	JKFF_Inst : JK_FlipFlop port map (J => '1', K => '1',
												 Clock => Q(i-1),
												 Enable => Enable,
												 Preset => Preset,
												 Reset => Reset,
												 Q => Q(i), NotQ => NotQ(i));
	end generate;
	
	Counter <= Q;
end Struct;