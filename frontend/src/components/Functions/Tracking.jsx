import React from 'react';
import ReactDOM from 'react-dom';
import { Redirect, withRouter } from "react-router-dom";
import { withStyles } from '@material-ui/core/styles';
import DropzoneComponent from 'react-dropzone-component';
import Button from '@material-ui/core/Button';
import CloudUploadIcon from '@material-ui/icons/CloudUpload';
import blue from '@material-ui/core/colors/blue';

import Paper from '@material-ui/core/Paper';

const styles = theme => ({
    buttoncss:{
        color: theme.palette.getContrastText(blue[600]),
        backgroundColor: blue[900],
        '&:hover': {
            backgroundColor:blue[600],
        },
        marginTop: '40px',
    },
    divCenter:{
        display:'flex',
        justifyContent:'center',
        alignItems:'center',
    },
})


class Tracking extends React.Component {

	constructor(props) {
		super(props);
		this.djsConfig = {
            dictDefaultMessage:"Drag a cgMLST profile here",
            dictRemoveFile:"Remove",
            addRemoveLinks: true,
            acceptedFiles: ".tsv",
            autoProcessQueue: false,
            parallelUploads: 1,
            maxFiles: 1,
            timeout: 0,
            init:function(){
                this.on("addedfile", function(file){
                    if(file.size < 1000){
                        this.removeFile(file);
                    }
                });
                this.on("sending", function(file, xhr, formData){
                    formData.append("profile_db", window.database);
                });
                this.on("success", function(file,response){
                    file._removeLink.remove();
                    delete file._removeLink;
                    window.trackingID = response.id;
                });
                this.on("maxfilesexceeded", function(file){
                    this.removeAllFiles();
                    this.addFile(file);
                });
            }
        }

        this.componentConfig = {
            iconFiletypes: ['.tsv'],
            showFiletypeIcon: true,
            postUrl: 'api/tracking/profile/'
        };

        this.dropzone = null;
	}

    handlePost() {
	    let fileCheck = this.dropzone.files.length;
	    if(fileCheck != 1){
	        alert('Please upload only 1 file');
	        return ;
	    }
	    this.dropzone.processQueue();
        this.props.history.push("/cgMLST/trackingResult")
	}

	render() {
		const config = this.componentConfig;
		const djsConfig = this.djsConfig;
        const eventHandlers = {
            init: dz => this.dropzone = dz,}
        const { classes } = this.props;

        return (
            <div style={{ marginTop:'100px' }}>
                { window.database == 'Salmonella_enterica' ? <Redirect to='/cgMLST/' /> : null }
                <div className={classes.divCenter}>
                    <DropzoneComponent config={config} eventHandlers={eventHandlers}
                        djsConfig={djsConfig} />
                </div>
                <div className={classes.divCenter}>
                    <Button variant="contained" className ={classes.buttoncss}
                     onClick={this.handlePost.bind(this)}>
                        Submit
                        &nbsp;&nbsp;
                        <CloudUploadIcon />
                    </Button>
                </div>
            </div>
        );
    }
}

export default withRouter(withStyles(styles)(Tracking));