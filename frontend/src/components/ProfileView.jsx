import React from 'react';
import ReactDOM from 'react-dom';
import { Link } from 'react-router-dom';
import Button from '@material-ui/core/Button';
import { withStyle } from '@material-ui/core/styles';
import ReplyIcon from '@material-ui/icons/Reply';
import DownloadIcon from '@material-ui/icons/CloudDownload';
import CircularProgress from '@material-ui/core/CircularProgress';
import download from 'downloadjs';

export default class Profile_view extends React.Component {

	constructor(props) {
		super(props);

        let fileName_str="", i=0, j=0;
        for (i; i < window.fileName.length; i++){
            if(i < window.fileName.length-1){
                fileName_str += window.fileName[i] + ", "
            }else{
                fileName_str += window.fileName[i]
            }
        };

        this.state = { fileName :fileName_str };
		this.query_result = this.query_result.bind(this);
    }

	query_result(){

		if(this.state.profile_result_zip == undefined){
			fetch('api/profiling/profile/' + window.batchid, { method:'GET'})
			.then(response => response.json())
			.then(result => this.setState(state => ({
                profile_result_zip: result.zip })));
            // fetch('api/dendrogram/dendrogram/' + window.batchid, { method: 'GET'})
            // .then(response => response.json())
            // .then(result => this.setState(state => ({
            //     png_file: result.png_file, 
            //     pdf_file: result.pdf_file,
            //     svg_file: result.svg_file, 
            //     emf_file: result.emf_file, 
            //     newick_file: result.newick_file })));
		}else{
			clearInterval(this.interval);
		}

	}

    getIdFile(){
        download(window.batchid,'BatchId.txt',"text/tab-separated-values");
    }

	componentDidMount(){
		this.query_result();
		this.interval = setInterval(this.query_result, 60000);
	}

    turn_on_Tabs(){
        window.tabSwitch = false;
    }

    render() {

    	if(this.state.profile_result_zip == undefined){

    		return(
    			<div>
                    <br />
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <font size="6"> ID : {window.batchid}</font>
                        &nbsp;&nbsp;&nbsp;&nbsp;
                        <Button variant="contained" color="default" onClick={this.getIdFile}>
                            Get ID
                        </Button>
                    </div>
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <font size="4">Database : {window.databaseName}</font>
                    </div>
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <font size="4">File name : {this.state.fileName}</font>
                    </div>
                    <br />
                    <br />
                    <br />
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <CircularProgress size={175} />
                    </div>
                    <br />
                    <br />
                    <br />
    			</div>
    			);
    	
    	}else{
    		return (
    				<div>
<<<<<<< HEAD
                        <br />
                        <br />
                        <br />
                        <br />
                        <br />
                        <br />
                        <br />
=======
                        <br />
                        <br />
    					<br />
                        <br />
                        <br />
                        <br />
                        <br />
                        <br />
                        <br />
                        <br />
                        <br />
                        <br />
>>>>>>> Modify user interface
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <a download href={this.state.profile_result_zip} 
                             style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Download profiles (.zip)
                                &nbsp;&nbsp;
                                <DownloadIcon />
                                </Button>
                            </a>
<<<<<<< HEAD
=======
                            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.profile_result_all} 
                             style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Download profiles (.tsv)
                                &nbsp;&nbsp;
                                <DownloadIcon />
                                </Button>
                            </a>
>>>>>>> Modify user interface
                        </div>
    					<br />
    					<br />
                        <br />
                        <br />
                        <br />
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <Link to="/" style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default" onClick={this.turn_on_Tabs}>
                                    <ReplyIcon />
                                    &nbsp;&nbsp;
                                    Back
                                </Button>
                            </Link>
                        </div>
                        <br />
                        <br />
    				</div>
        	);
    	}
        
    }
<<<<<<< HEAD
}
=======
}

                        // <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        //     <img src={this.state.svg_file} />
                        // </div>
                        // <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        //     <font>Download</font> 
                        //     &nbsp;&nbsp;&nbsp;&nbsp;
                        //     <a download href={this.state.pdf_file} style={{ textDecoration:'none' }}>
                        //         <Button variant="contained" color="default">
                        //         Pdf
                        //         </Button>
                        //     </a>
                        //     &nbsp;&nbsp;&nbsp;&nbsp;
                        //     <a download href={this.state.emf_file} style={{ textDecoration:'none' }}>
                        //         <Button variant="contained" color="default">
                        //         emf
                        //         </Button>
                        //     </a>
                        //     &nbsp;&nbsp;&nbsp;&nbsp;
                        //     <a download href={this.state.svg_file} style={{ textDecoration:'none' }}>
                        //         <Button variant="contained" color="default">
                        //         Svg
                        //         </Button>
                        //     </a>
                        //     &nbsp;&nbsp;&nbsp;&nbsp;
                        //     <a download href={this.state.png_file} style={{ textDecoration:'none' }}>
                        //         <Button variant="contained" color="default">
                        //         Png 
                        //         </Button>
                        //     </a>
                        //     &nbsp;&nbsp;&nbsp;&nbsp;
                        //     <a download href={this.state.newick_file} style={{ textDecoration:'none' }}>
                        //         <Button variant="contained" color="default">
                        //         newick
                        //         </Button>
                        //     </a>
                        // </div>
>>>>>>> Modify user interface
