import React from 'react';
import ReactDOM from 'react-dom';
import { withStyles } from '@material-ui/core/styles';
import { Link } from 'react-router-dom';
import { Prompt } from 'react-router';
import Button from '@material-ui/core/Button';
import ReplyIcon from '@material-ui/icons/Reply';
import CircularProgress from '@material-ui/core/CircularProgress';
import download from 'downloadjs';

const styles = theme => ({
    displayCenter:{
        display:'flex',
        justifyContent:'center',
        alignItems:'center',
    },
    button:{
        marginLeft: '25px',
    },
})

class ClusterinhResult extends React.Component {

	constructor(props) {
		super(props);
        this.state = {};
		this.queryDengrogram = this.queryDengrogram.bind(this);
	}

    componentDidMount(){
		this.interval = setInterval(this.queryDengrogram, 1000);
	}

    downloadJobID(){
        download(window.clusteringID,'clustering_JobID.txt',"text/tab-separated-values");
    }

    queryDengrogram(){

		if(this.state.pngFile == undefined){
			fetch('api/dendrogram/dendrogram/' + window.clusteringID, { method: 'GET'})
			.then(response => response.json())
			.then(result => this.setState(state => ({
				pngFile: result.png_file,
				pdfFile: result.pdf_file,
				svgFile: result.svg_file,
				newickFile: result.newick_file, })));

		}else {
			clearInterval(this.interval);
		}

	}

	render() {
		const { classes } = this.props;

    	if(this.state.pngFile == undefined){
    		return(
	    		<div style={{ marginTop: '100px', }}>
	    			<Prompt
                            when={true}
                            message="You are leaving the page. Please save ID to get result. Are you sure to leave now?"/>
					<div className={classes.displayCenter}>
                        <font size="6"> Job ID : {window.clusteringID}</font>
                        &nbsp;&nbsp;&nbsp;&nbsp;
                        <Button variant="contained" color="default" onClick={this.downloadJobID}>
                            Get job ID
                        </Button>
                    </div>
					<div style={{ marginTop: '20px', }} className={classes.displayCenter}>
						<font size="6"> Please hold on ... </font>
					</div>

					<div style={{ marginTop: '100px', marginBottom:'100px' }} className={classes.displayCenter}>
						<CircularProgress size={175} />
	                </div>
				</div>
		    );
    	}else{
    		return (
				<div style={{ marginTop: '100px', }}>
					<Prompt
                            when={true}
                            message="You are leaving the page. Please save results, or it will lose. Are you sure to leave now?"/>
					<div className={classes.displayCenter}>
                        <font size="6">Job ID : {window.clusteringID}</font>
                        <Button style={{ marginLeft: '20px', }} variant="contained" color="default" onClick={this.downloadJobID}>
                            Get job ID
                        </Button>
                    </div>
					<div style={{ margin: '20px', }} className={classes.displayCenter}>
						<img src={this.state.svgFile} />
					</div>
					<br />
					<div className={classes.displayCenter}>
						<font>Download</font>
					</div>
					<br />
					<div className={classes.displayCenter}>
						<a download href={this.state.pngFile} style={{ textDecoration:'none' }}>
							<Button className={classes.button} variant="contained" color="default">Png</Button>
						</a>
						<a download href={this.state.pdfFile} style={{ textDecoration:'none' }}>
							<Button className={classes.button} variant="contained" color="default">Pdf</Button>
						</a>
						<a download href={this.state.svgFile} style={{ textDecoration:'none' }}>
							<Button className={classes.button} variant="contained" color="default">Svg</Button>
						</a>
						<a download href={this.state.newickFile} style={{ textDecoration:'none' }}>
							<Button className={classes.button} variant="contained" color="default">newick</Button>
						</a>
					</div>
					<div className={classes.displayCenter}>
						<Link to="/cgMLST/clustering" style={{ textDecoration:'none' }}>
							<Button style={{ marginTop: '50px' }} variant="contained" color="default">
								<ReplyIcon />
								&nbsp;&nbsp;
								Back
							</Button>
						</Link>
					</div>
				</div>
			);
		}
    }
}

export default withStyles(styles)(ClusterinhResult);