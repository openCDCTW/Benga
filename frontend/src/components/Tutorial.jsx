import React from 'react';
import ReactDOM from 'react-dom';
import { Link } from 'react-router-dom';
import { withStyles } from '@material-ui/core/styles';
import Paper from '@material-ui/core/Paper';
import Typography from '@material-ui/core/Typography';

const styles = theme => ({
	paper_root: {
		...theme.mixins.gutters(),
		paddingTop: theme.spacing.unit * 2,
		paddingBottom: theme.spacing.unit * 3,
		width:'90%',
	},
})

class Tutorial extends React.Component {

	render() {
		const { classes } = this.props;

		return (
			<div>
				<br />
				<div style={{ justifyContent:'center',margin:'auto',display:'flex'}}>
					<Paper className={classes.paper_root} elevation={4}>
						<h2>&nbsp;How to use?</h2>
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto',fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;When connect to this page, 
							you will see the following interface.
						</Typography >
						<br />
						<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
							<img style={{ width:"90%" }} src={require('./static/default_page.png')}/>
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto',fontSize:16}}>
						&nbsp;&nbsp;&nbsp;&nbsp;Different service from navigation bar.
						</Typography >
						<br />
						<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
							<img style={{ width:"90%" }} src={require('./static/navigation.gif')} />
						</div>
					</Paper>
				</div>
				<br />
				<div style={{ justifyContent:'center',margin:'auto',display:'flex'}}>
					<Paper className={classes.paper_root} elevation={4}>
						<h2>&nbsp;cgMLST Profiling</h2>
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto',fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;You can upload sequences to get cgMLST profiles 
							at this pagination. You can drop files (or cilck) to specific area to add 
							files to be upload as shown in figure. Only accept file format of .fa, .fna, 
							.fasta.
						</Typography >
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img style={{ width:"90%" }} src={require('./static/dropzone.gif')} />
						</div>					
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto',fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;When you drag file(s) to dropzone, 
							it will display file name, file size, 
							progress bar and remove link as shown in fingure.
						</Typography >
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img style={{ width:"60%" }} src={require('./static/fileInfo.png')} />
						</div>	
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto',fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;If you wants to remove specified file, click the 
							“remove” link. If you wants to remove all files, click “remove all files” button.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
							<img style={{ width:"90%" }} src={require('./static/remove.gif')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto',fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;After drag file(s), choose a cgMLST database.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img style={{ width:"40%" }} src={require('./static/select.gif')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp; Then click the "submit" button to upload file(s) and do profiling.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img style={{ width:"90%" }} src={require('./static/profile_submit.gif')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;After submit, the screen will display job ID, database and file name(s).
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img style={{ width:"70%" }} src={require('./static/profiling_info.png')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;Now you can wait for completion, or copy(download) job ID 
							result. Result page is shown as figure.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img style={{ width:"90%" }} src={require('./static/profile_result.png')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;After process complete, you can input job ID in seach bar. 
							Then press enter, the screen will show the result.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img style={{ width:'90%'}} src={require('./static/profile_ID_result.gif')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;If the process still ongoing or input wrong job ID, 
							the warning will shown.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img style={{ width:'90%'}} src={require('./static/profile_ID_result_fail.gif')} />
						</div>
						<br />
					</Paper>
				</div>
				<br />
				<div style={{ justifyContent:'center',margin:'auto',display:'flex'}}>
					<Paper className={classes.paper_root} elevation={4}>
						<h2>&nbsp;&nbsp;Strain Tracking</h2>
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;You can upload a cgMLST profile to compare cgMLST 
							profiles of global outbreak strains. We only provide Vibrio Cholerae database in this vision. 
							Before tracking, you have to upload a cgMLST profile (only accept .tsv file format). 
							After upload file, click “submit” button to do tracking. 
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img style={{ width:'90%' }} src={require('./static/tracking_submit.gif')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;Then wait for completion. The result will shown as figure.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img style={{ width:'90%' }} src={require('./static/tracking_done.gif')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;You can click "download" button to get table and profiles. 
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img style={{ width:'90%' }} src={require('./static/tracking_result.png')} />
						</div>
						<br />
						<br />
					</Paper>
				</div>
				<br />
				<div style={{ justifyContent:'center',margin:'auto',display:'flex'}}>
					<Paper className={classes.paper_root} elevation={4}>
						<h2>&nbsp;&nbsp;Clustering</h2>
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;You can upload profiles to get 
							dendrogram. You have to upload at least 5 cgMLST profiles(only accept .tsv).
						</Typography>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;You can select singale linkage or UPGMA algorithm to 
							draw dendrogram.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img style={{ width:'50%' }} src={require('./static/clustering_select.gif')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;After upload files, click “submit” button to draw dendrogram.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img style={{ width:'90%' }} src={require('./static/clustering_submit.gif')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							After completion, you can download dendrogram of png, pdf, svg, emf and newick format. 
							We also provide using job ID to get result. The job ID is valid for two weeks.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img style={{ width:'90%' }} src={require('./static/clustering_result.png')} />
						</div>
						<br />
					</Paper>
				</div>
				<br />
			</div>
		)
	}

}

export default withStyles(styles)(Tutorial)