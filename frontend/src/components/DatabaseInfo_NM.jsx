import React from 'react';
import ReactDOM from 'react-dom';
import Typography from '@material-ui/core/Typography';
import { withStyles } from '@material-ui/core/styles';
import Button from '@material-ui/core/Button';
import DownloadIcon from '@material-ui/icons/CloudDownload';
import green from '@material-ui/core/colors/green';
import NM_img from './static/Nm_loci_feq_distribution.png';

const styles = theme => ({
        downloadButton:{
          color: theme.palette.getContrastText(green[900]),
          backgroundColor: green[600],
          '&:hover': {
            backgroundColor:green[900],
          },
        },
});

class DatabaseInfo_NM extends React.Component {

	constructor(props) {
		super(props);
	}

	render(){
		const { classes } = this.props;

		return (
			<div>
				<Typography component="p" style={{ width:'90%',textAlign: 'justify',
				margin:'auto',fontSize:16 }}>
					cgMLST@Taiwan provides a 
					<font style={{ fontStyle:'italic' }}> Neisseria meningitidis </font> 
					allele database,
					<font style={{ fontStyle:'italic' }}> N. meningitidis </font>
					cgMLST profile database, and tools for cgMLST profiling, strain tracking, 
					and clustering of cgMLST profiles via the internet. cgMLST profiling 
					is based on 1,241 
					<font style={{ fontStyle:'italic' }}> N. meningitidis </font> 
					core genes, which are identified from 319
					<font style={{ fontStyle:'italic' }}> N. meningitidis </font> 
					genomes obtained from the NCBI database. Core genes are designated 
					for those existing in more than 95% of the 319 genomes.
				</Typography>
				<br />
				<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                	<img style={{ width:'95%' }} src={NM_img} />
            	</div>
            	<br />
            	<Typography component="p" style={{ width:'90%',
				margin:'auto',fontSize:18, fontWeight:500 }}>
					Figure. Frequency of loci (genes) over 319
					<font style={{ fontStyle:'italic' }}> N. meningitidis </font>
					genomes.
				</Typography>
				<br />
                <br />
				<Typography component="p" style={{ width:'90%',textAlign: 'justify',
				margin:'auto',fontSize:16 }}>
					The cgMLST profile database contains cgMLST profiles for the 
					<font style={{ fontStyle:'italic' }}> N. meningitidis </font> 
					strains with genomic sequences deposited in the NCBI database. 
					Nowadays, the database contains 9,431 cgMLST profiles and will be updated 
					over time.
				</Typography>
				<br />
				<Typography component="p" style={{ width:'90%',textAlign: 'justify',
				margin:'auto',fontSize:16 }}>
					cgMLST@Taiwan was developed in the laboratories of Centers for Disease Control, 
					Ministry of Health and Welfare, Taiwan by Yueh-Hua Tu, Yi-Syong Chen, Bo-Han Chen, 
					Yen-Yi Liu, Yu-Ping Hong, Ru-Hsiou Teng, You-Wun Wang, and Chien-Shun Chiou.
				</Typography>
                <br />
                <div style = {{ width:'90%', margin:'auto' }}>
                  <a download href='https://drive.google.com/uc?export=download&id=1sYZ-dnUZwYLG4ZgjXDPiLDBdgHiM5Zvh' style={{ textDecoration:'none' }}>
                     <Button style={{ textTransform:'none' }} variant="contained" className={classes.downloadButton}>
                        Download &nbsp; 1,241 &nbsp; N. meningitidis &nbsp; core &nbsp; genes &nbsp; (.csv) &nbsp;
                        <DownloadIcon />
                     </Button>
                  </a>
                </div>
				<br />
				<h3 style={{ marginLeft:'30px' }}>Contact</h3>
				<Typography component="p" style={{ width:'90%',
				margin:'auto',fontSize:16 }}>
					For any question please contact Chien-Shun Chiou by email:
				</Typography>
				<Typography component="p" style={{ width:'90%',textAlign: 'justify',
				margin:'auto',fontSize:16 }}>
					<font style={{ color:'blue', textDecoration:'underline' }}>nipmcsc@cdc.gov.tw</font>
					 &nbsp; or &nbsp;
					 <font style={{ color:'blue', textDecoration:'underline' }}>nipmcsc@gmail.com</font>
				</Typography>
				<br />
			</div>
		)
	}

}
export default withStyles(styles)(DatabaseInfo_NM);