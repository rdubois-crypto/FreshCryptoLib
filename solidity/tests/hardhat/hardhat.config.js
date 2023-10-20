require("@nomicfoundation/hardhat-toolbox");
require("hardhat-gas-reporter");
require("hardhat-preprocessor");
const removeConsoleLog = require("hardhat-preprocessor").removeConsoleLog;

const fs = require("fs");
const mnemonic = fs.readFileSync(".mnemonic").toString().trim()

task("accounts", "Prints the list of accounts", async (taskArgs, hre) => {
  const accounts = await hre.ethers.getSigners();

  for (const account of accounts) {
    console.log(account.address);
  }
});

/** @type import('hardhat/config').HardhatUserConfig */
module.exports = {
  solidity: "0.8.17",
   settings: {
       viaIR: true, // add this

      optimizer: {
        enabled: true,
        runs: 100000,
      },
    },
    
	preprocess: {
    eachLine: removeConsoleLog((hre) => hre.network.name !== 'hardhat' && hre.network.name !== 'localhost'),
  },    
  defaultNetwork: "hardhat",
  gasReporter: {
    enabled: (process.env.REPORT_GAS) ? true : false
  },
  networks: {
  	hardhat: {
/*      
			forking: { 
				url: "http://192.168.1.4:8545",
				blockNumber: 16665200
			}
*/            
  	},
  	polygon_zk_test: {
  		url : "https://rpc.public.zkevm-test.net",
  		accounts: {
  			mnemonic: mnemonic  			
  		}
  	},
  	goerli: {
  		url : "https://endpoints.omniatech.io/v1/eth/goerli/public",
  		accounts: {
  			mnemonic: mnemonic  			
  		}  		
  	}
  }
};
