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
import ReplyIcon from '@material-ui/icons/Reply';
import DownloadIcon from '@material-ui/icons/CloudDownload';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import blue from '@material-ui/core/colors/blue';
import download from 'downloadjs';


const styles = theme => ({
    cssRoot:{
        color: theme.palette.getContrastText(blue[900]),
        backgroundColor: blue[900],
        '&:hover': {
            backgroundColor:blue[600],
        },
    },
    divCenter:{
        display:'flex',
        justifyContent:'center',
        alignItems:'center',
    },
})

class Profiling extends React.Component {

    constructor(props) {
        super(props);
        window.fileName = [];

        this.djsConfig = {
            dictDefaultMessage:"Drag contig file(s) here (up to 100 files)",
            dictRemoveFile:"Remove",
            addRemoveLinks: true,
            maxFilesize: 10,
            acceptedFiles: ".fasta,.fa,.fna",
            autoProcessQueue: false,
            timeout: 0,
            parallelUploads: 100,
            init:function(){
                this.on("addedfile", function(file){
                    if(file.size < 10000){
                        this.removeFile(file);
                    }
                });
                this.on("sending", function(file, xhr, formData){
                    formData.append("batch_id", window.batchid);
                    formData.append("database", window.database);
                    formData.append("occurrence", window.occurrence);
                    window.fileName.push(file.name);
                });
                this.on("success", function(file){
                    file._removeLink.remove();
                    delete file._removeLink;
                });
            }
        }

        this.componentConfig = {
            iconFiletypes: ['.fasta','.fna','.fa'],
            showFiletypeIcon: true,
            postUrl: 'api/profiling/sequence/'
        };

        this.dropzone = null;
    }

    handlePost() {

        var fileCheck = this.dropzone.files.length;

        if(fileCheck < 1){
            alert('Please upload contigs file first');
            return ;
        }else if(fileCheck > 100){
            alert('Cannot upload more than 100 files');
            return ;
        }

        fetch('api/profiling/upload/', {method:'POST'})
            .then(function(res){
               return res.json();
            }).then(batch => window.batchid = batch.id);

        function submit(){

            if( window.batchid != undefined ){
                this.dropzone.processQueue();

                var scheme = {};
                scheme.seq_num = fileCheck;
                fetch('api/profiling/upload/' + window.batchid + '/', {
                    method:'PATCH',
                    headers: new Headers({'content-type': 'application/json'}),
                    body: JSON.stringify(scheme)
                });

                this.props.history.push("/cgMLST/profilingResult");
                clearInterval(interval);
            }
        };

        let interval = setInterval(submit.bind(this),500);
    }

    remove(){
        this.dropzone.removeAllFiles();
        window.fileName.length = 0;
    }

    query(){

        fetch('api/profiling/profile/' + this.state.querybyID + '/', { method:'GET'})
        .then(function(response){
            if(response.status != 404){
                return response.json();
            }else{
                return response.status;
            }
        }).then(res => this.setState(state => ({ tmp: res })));

        function result(){
            if(this.state.tmp.zip != undefined){
                this.setState(state => ({ profile_result_zip : this.state.tmp.zip }))
            }else{
                alert("Data not found. Please input correct ID or try again later.");
            }
            clearInterval(interval);
        }
        let interval = setInterval(result.bind(this),50);
    }

    render() {

        const config = this.componentConfig;
        const djsConfig = this.djsConfig;
        const eventHandlers = {
            init: dz => this.dropzone = dz,}
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
                <br />
                <div className={classes.divCenter}>
                    { window.database == 'Vibrio_cholerae' ?
                    <a download href='https://drive.google.com/uc?export=download&id=1XG-05Kim8gOOg1UU16oJ9saOwfCh-PZM'
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

export default withRouter(withStyles(styles)(Profiling));