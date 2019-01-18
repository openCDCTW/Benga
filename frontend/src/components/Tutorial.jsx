import React from 'react';
import ReactDOM from 'react-dom';
import Navigation from './Navigation.jsx';
import { Link } from 'react-router-dom';


export default class Tutorial extends React.Component {

	render() {
		return (
			<div>
				<Navigation value={4}/>
				<br />
				<br />
				<h2>&nbsp;&nbsp; About this website</h2>
				<div style={{ width:'90%'}}>
					<font>
						&nbsp;&nbsp;&nbsp;&nbsp;It’s a web-basde analysis tool for bacterial 
						genomes. This website has three main features. One is upload sequences 
						to get cgMLST profiles as well as dendrogram. Another is upload profile 
						to get dendrogram. The other is upload sequences to compare cgMLST 
						profiles of global outbreak strains. 
					</font>
				</div>
				<br/>
				<h2>&nbsp;&nbsp; How to use?</h2>
				<div style={{ width:'90%'}}>
					<font>
						&nbsp;&nbsp;&nbsp;&nbsp;When user connect to this page, user will see 
						the following interface if there is no accident.
					</font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<img src={require('./static/default_page.png')}/>
				</div>
				<br />
				<div style={{ width:'90%'}}>
					<font>
						&nbsp;&nbsp;&nbsp;&nbsp;User can choose different service from navigation bar.
					</font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<img src={require('./static/navigation.png')} />
				</div>
				<br />
				<hr color = "#000000" size="5" />
				<h3>&nbsp;&nbsp;Profiling</h3>
				<br />
				<div style={{ width:'90%' }}>
					<font>
						&nbsp;&nbsp;&nbsp;&nbsp;At profiling pagination, user can upload sequences to get cgMLST 
						profiles as well as dendrogram. User can drop files (or cilck) to specific 
						area to add files to be upload as shown in figure. Before doing profiling 
						and draw dendrogram, user have to upload at least 5 sequence files 
						(only accept .fa, .fna, .fasta).
					</font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<img src={require('./static/dropzone.png')} />
				</div>
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<font>
						Before add files
					</font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<img src={require('./static/dropzone_addedfiles.png')} />
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<font>
						After add files
					</font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<font>
						&nbsp;&nbsp;&nbsp;&nbsp;If user wants to remove specified file, user can click “remove file” 
						link below the thumbnail.
					</font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<img src={require('./static/Thumbnail.png')} />
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center' }}>
					<font>
						&nbsp;&nbsp;&nbsp;&nbsp;If user wants to remove all file, user can click “remove all files” 
						button on the screen.					
					</font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<img src={require('./static/remove_button.png')} />
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<font>
						&nbsp;&nbsp;&nbsp;&nbsp;After drop sequences files to specific area, choose a cgMLST 
						database.
					</font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<img src={require('./static/select.png')} />
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<font>
						&nbsp;&nbsp;&nbsp;&nbsp;Then click the upload button to upload files. Now, user can click 
						profiling button to do profiling and draw dendrogram.
					</font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<img src={require('./static/upload_profiling_button.png')} />
				</div>
				<br />
				<div style={{ width:'90%' }}>
					<font>
						&nbsp;&nbsp;&nbsp;&nbsp;TAfter click “PROFILING” button, the screen will display user’s 
						upload information as shown in figure.
					</font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<img src={require('./static/doing_profiling.png')} />
				</div>
				<br />
				<div style={{ width:'90%' }}>
					<font>
						&nbsp;&nbsp;&nbsp;&nbsp;Users can choose not to close the page and wait for completion. 
						Or user can use batch ID to query result if the process has been done.User 
						can input batch ID in the search bar. Then press enter, the screen will 
						show the result.
					</font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<img src={require('./static/seachbar.png')} />
				</div>
				<br />
				<hr color = "#000000" size="5" /> 
				<h3>&nbsp;&nbsp;Dendrogram</h3>
				<div style={{ width:'90%' }}>
					<font>
						&nbsp;&nbsp;&nbsp;&nbsp;At dendrogram pagination, user can upload profile to get 
						dendrogram. Upload steps are as same as profiling pagination. Before 
						draw dendrogram, user have to upload at least 5 profile files 
						(only accept .tsv).  After upload files, user can click “draw dendrogram” 
						button to draw dendrogram. 
					</font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<img src={require('./static/upload_dendrogram_button.png')} />
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<font>
						Then wait for completion after click “DRAW DENDROGRAM” button.
					</font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<img src={require('./static/waiting.png')} />
				</div>
				<br />
				<hr color = "#000000" size="5" /> 
				<h3>&nbsp;&nbsp;Tracking</h3>
				<div style={{ width:'90%'}}>
					<font>
						&nbsp;&nbsp;&nbsp;&nbsp;At dendrogram pagination, upload sequences to compare cgMLST 
						profiles of global outbreak strains. Upload steps are as same as before. 
						We only provide Vibrio Cholerae database in this vision. Before tracking, 
						user have to upload at least 1 sequence file (only accept .fa, .fna, .fasta).  After upload files, user can click “tracking” button to do tracking.
					</font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<img src={require('./static/upload_tracking_button.png')} />
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<font>
						Then wait for completion after click “TRACKING” button. 
					</font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', 
				alignItems:'center'}}>
					<img src={require('./static/waiting.png')} />
				</div>
				<br />
				<br />
			</div>
		)
	}

}
