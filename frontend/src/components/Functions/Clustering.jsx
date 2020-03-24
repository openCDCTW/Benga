import React from 'react';
import ReactDOM from 'react-dom';
import { withRouter } from "react-router-dom";
import DropzoneComponent from 'react-dropzone-component';
import { Link } from 'react-router-dom';
import Button from '@material-ui/core/Button';
import { withStyles } from '@material-ui/core/styles';
import CloudUploadIcon from '@material-ui/icons/CloudUpload';
import Icon from '@material-ui/core/Icon';
import DeleteIcon from '@material-ui/icons/Delete';
import blue from '@material-ui/core/colors/blue';
import Radio from '@material-ui/core/Radio';
import ReplyIcon from '@material-ui/icons/Reply';

const styles = theme => ({
    cssRoot:{
        color: theme.palette.getContrastText(blue[600]),
        backgroundColor: blue[900],
        '&:hover': {
            backgroundColor:blue[600],
        },
    },
    radio:{
        '&$checked': {
            color: blue[600],
        },
    },
    divCenter:{
        display:'flex',
        justifyContent:'center',
        alignItems:'center',
    },
})

class Clustering extends React.Component {

    constructor(props) {
        super(props);
        this.state = {
            algorithm: 'single'
        };

        this.djsConfig = {
            dictDefaultMessage:"Drag cgMLST profiles here (up to 500 files)",
            dictRemoveFile:"Remove",
            addRemoveLinks: true,
            acceptedFiles: ".tsv",
            autoProcessQueue: false,
            parallelUploads: 500,
            timeout: 0,
            init:function(){
                this.on("addedfile", function(file){
                    if(file.size < 1000){
                        this.removeFile(file);
                    }
                });
                this.on("sending", function(file, xhr, formData){
                    formData.append("batch_id", window.clusteringID);
                });
                this.on("success", function(file){
                    file._removeLink.remove();
                    delete file._removeLink;
                });
            }
        }

        this.componentConfig = {
            iconFiletypes: ['.tsv'],
            showFiletypeIcon: true,
            postUrl: 'api/dendrogram/profile/'
        };

        this.dropzone = null;
    }

    handlePost() {
        let fileCheck = this.dropzone.files.length;
        if(fileCheck < 1){
            alert('Please upload at least 1 files');
            return ;
        }else if(fileCheck > 500){
            alert('Cannot upload more than 500 files');
            return ;
        }

        fetch('api/dendrogram/upload/', {method:'POST'})
            .then(function(res){
               return res.json();
            }).then(batch => window.clusteringID = batch.id);

        function submit(){
            if( window.clusteringID != undefined ){
                this.dropzone.processQueue();
                let scheme = {};
                scheme.prof_num = this.dropzone.files.length;
                scheme.linkage = this.state.algorithm;
                fetch('api/dendrogram/upload/' + window.clusteringID + '/', {
                    method:'PATCH',
                    headers: new Headers({'content-type': 'application/json'}),
                    body: JSON.stringify(scheme)
                });
                this.props.history.push("/cgMLST/clusteringResult");
                clearInterval(interval);
            }
        };
        let interval = setInterval(submit.bind(this),500);
    }

    remove(){
        this.dropzone.removeAllFiles();
    }

    handleChange(event){
        this.setState({ algorithm: event.target.value });
    }

    render() {
        const config = this.componentConfig;
        const djsConfig = this.djsConfig;
        const eventHandlers = {
            init: dz => this.dropzone = dz,
        }
        const { classes } = this.props;

        return (
            <div style={{ marginTop:'100px' }}>
                <div style={{ float:'right', marginTop:'-50px', marginRight:'25px' }}>
                    <Button variant="contained" color="secondary" onClick={this.remove.bind(this)}>
                        Remove all files
                        &nbsp;&nbsp;
                        <DeleteIcon />
                    </Button>
                </div>
                <div className={classes.divCenter}>
                    <DropzoneComponent config={config} eventHandlers={eventHandlers}
                    djsConfig={djsConfig} />
                </div>
                <br />
                <div className={classes.divCenter}>
                    <font>Algorithms : </font>
                    &nbsp;&nbsp;
                    <Radio
                        color='primary'
                        checked={this.state.algorithm === 'single'}
                        onChange={this.handleChange.bind(this)}
                        value="single"
                    />
                    <font>Single linkage</font>
                    <Radio
                        color='primary'
                        checked={this.state.algorithm === 'average'}
                        onChange={this.handleChange.bind(this)}
                        value="average"
                    />
                    <font>UPGMA</font>
                </div>
                <br />
                <br />
                <div className={classes.divCenter}>
                    { window.database == 'Vibrio_cholerae' ?
                        <a download href='https://drive.google.com/uc?export=download&id=1RyJl1yMrWRHPgYAbStjUtD2zmfQ8YMRb'
                        style={{ textDecoration:'none', marginRight:'25px' }}>
                            <Button style={{ textTransform:'none' }} variant="contained" color="default">
                                Download &nbsp; example &nbsp; files
                            </Button>
                        </a> : null }
                    <Button variant="contained" className ={classes.cssRoot}
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

export default withRouter(withStyles(styles)(Clustering));
