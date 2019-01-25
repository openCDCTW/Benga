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
						<h2>&nbsp;About this website</h2>			
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto',fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;It’s a web-basde analysis tool for bacterial 
							genomes. This website has three main features. One is upload sequences 
							to get cgMLST profiles as well as dendrogram. Another is upload profile 
							to get dendrogram. The other is upload sequences to compare cgMLST 
							profiles of global outbreak strains. 
						</Typography>
					</Paper>
				</div>
				<br/>
				<div style={{ justifyContent:'center',margin:'auto',display:'flex'}}>
					<Paper className={classes.paper_root} elevation={4}>
						<h2>&nbsp;How to use?</h2>
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto',fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;When user connect to this page, user will see 
							the following interface if there is no accident.
						</Typography >
						<br />
						<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
							<img src={require('./static/default_page.png')}/>
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto',fontSize:16}}>
						&nbsp;&nbsp;&nbsp;&nbsp;User can choose different service from navigation bar.
						</Typography >
						<br />
						<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
							<img src={require('./static/navigation.png')} />
						</div>
					</Paper>
				</div>
				<br />
				<div style={{ justifyContent:'center',margin:'auto',display:'flex'}}>
					<Paper className={classes.paper_root} elevation={4}>
						<h2>&nbsp;Profiling</h2>
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto',fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;At profiling pagination, user can upload sequences to get cgMLST 
							profiles as well as dendrogram. User can drop files (or cilck) to specific 
							area to add files to be upload as shown in figure. Before doing profiling 
							and draw dendrogram, user have to upload at least 5 sequence files 
							(only accept .fa, .fna, .fasta).
						</Typography >
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img src={require('./static/dropzone.png')} />
						</div>
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', display:'flex', justifyContent:'center',fontSize:16}}>
							Before add files
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img src={require('./static/dropzone_addedfiles.png')} />
						</div>
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto',display:'flex', justifyContent:'center',fontSize:16}}>
							After add files
						</Typography>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto',fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;If user wants to remove specified file, user can 
							click “remove file” link below the thumbnail.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
							<img src={require('./static/Thumbnail.png')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto',fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;If user wants to remove all file, 
							user can click “remove all files” button on the screen.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img src={require('./static/remove_button.png')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;After drop sequences files to specific area, 
							choose a cgMLST database.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img src={require('./static/select.png')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;Then click the upload button to upload files. 
							Now, user can click profiling button to do profiling and draw dendrogram.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img src={require('./static/upload_profiling_button.png')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;TAfter click “PROFILING” button, 
							the screen will display user’s upload information as shown in figure.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img src={require('./static/doing_profiling.png')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;Users can choose not to close the page and 
							wait for completion. Or user can use batch ID to query result if the 
							process has been done.User can input batch ID in the search bar. 
							Then press enter, the screen will show the result.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img style={{ width:'95%'}} src={require('./static/seachbar.png')} />
						</div>
						<br />
					</Paper>
				</div>
				<br />
				<div style={{ justifyContent:'center',margin:'auto',display:'flex'}}>
					<Paper className={classes.paper_root} elevation={4}>
						<h2>&nbsp;&nbsp;Dendrogram</h2>
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;At dendrogram pagination, user can upload profile to get 
							dendrogram. Upload steps are as same as profiling pagination. Before 
							draw dendrogram, user have to upload at least 5 profile files 
							(only accept .tsv).  After upload files, user can click “draw dendrogram” 
							button to draw dendrogram.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img src={require('./static/upload_dendrogram_button.png')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							Then wait for completion after click “DRAW DENDROGRAM” button.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img src={require('./static/waiting.png')} />
						</div>
						<br />
					</Paper>
				</div>
				<br />
				<div style={{ justifyContent:'center',margin:'auto',display:'flex'}}>
					<Paper className={classes.paper_root} elevation={4}>
						<h2>&nbsp;&nbsp;Tracking</h2>
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							&nbsp;&nbsp;&nbsp;&nbsp;At dendrogram pagination, upload sequences to compare cgMLST 
							profiles of global outbreak strains. Upload steps are as same as before. 
							We only provide Vibrio Cholerae database in this vision. Before tracking, 
							user have to upload at least 1 sequence file (only accept .fa, .fna, .fasta).  
							After upload files, user can click “tracking” button to do tracking.
						</Typography>
						<br />
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img src={require('./static/upload_tracking_button.png')} />
						</div>
						<br />
						<Typography component="p" style={{ width:'90%',textAlign: 'justify',
						margin:'auto', fontSize:16 }}>
							Then wait for completion after click “TRACKING” button.
						</Typography>
						<div style={{ display:'flex', justifyContent:'center', 
						alignItems:'center'}}>
							<img src={require('./static/waiting.png')} />
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