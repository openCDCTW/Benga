import React from 'react';
import ReactDOM from 'react-dom';
import DropzoneComponent from 'react-dropzone-component';
import { Link } from 'react-router-dom';
import Button from '@material-ui/core/Button';
import { withStyles } from '@material-ui/core/styles';
import CloudUploadIcon from '@material-ui/icons/CloudUpload';
import Icon from '@material-ui/core/Icon';
import DeleteIcon from '@material-ui/icons/Delete';
import green from '@material-ui/core/colors/green';
import blue from '@material-ui/core/colors/blue';
//
import Radio from '@material-ui/core/Radio';

const styles = theme => ({
    cssRoot:{
        color: theme.palette.getContrastText(green[600]),
        backgroundColor: green[500],
        '&:hover': {
            backgroundColor:green[600],
        },
    },
    radio:{
        '&$checked': {
            color: blue[600],
    },
},

})

class Upload_profile extends React.Component {

    constructor(props) {

        super(props);

        fetch('api/dendrogram/upload/', {method:'POST'})
        .then(function(res){
           return res.json();
        }).then(function(batch){
           return getID(batch);
        });

        var getID=function(data){
            window.clusteringID = data.id;
        };

        // TODO: poor performnce issue with everytime load the component will
        // fetch API once.

        this.state = {
            to: "/upload_profile",
            upload_confirm: false,
            algorithm: 'singal_linkage'
        };

        this.djsConfig = {
            dictDefaultMessage:"Drop or click to upload profiles",
            addRemoveLinks: true,
            acceptedFiles: ".tsv",
            autoProcessQueue: false,
            parallelUploads: 500,
            init:function(){
                this.on("addedfile", function(on_load_header_data){

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

        this.dropzone.processQueue();
        this.setState(state => ({ to: '/dendrogram_view', upload_confirm: true}));
    
    }

    submit(){

        if (this.state.upload_confirm == false ){
            alert('Please upload files first! (At least 1 files)');
            return ;
        }

        var scheme = {};
        scheme.id = window.clusteringID;
        scheme.linkage = this.state.algorithm;
        fetch('api/dendrogram/plot/', {
            method:'POST',
            headers: new Headers({'content-type': 'application/json'}),
            body: JSON.stringify(scheme)
        });

        window.tabSwitch = true;
    }

    remove(){

        this.dropzone.removeAllFiles();
        
        this.setState(state => ({ to: '/upload_profile', upload_confirm: false}));
        fetch('api/dendrogram/upload/', {method:'POST'})
        .then(function(res){
           return res.json();
        }).then(function(batch){
           return getID(batch);
        });

        var getID=function(data){
            window.clusteringID = data.id;
        };

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
            <div>
                <div>
                    <div style={{ float:'right', marginTop:'35px', marginRight:'25px' }}>
                        <Button variant="contained" color="secondary" onClick={this.remove.bind(this)}>
                            Remove all files
                            &nbsp;&nbsp;
                            <DeleteIcon />
                        </Button>
                    </div>
                </div>
                <br />
                <br />
                <br />
                <br />
                <br />
                <div style = {{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                    <DropzoneComponent config={config} eventHandlers={eventHandlers} 
                    djsConfig={djsConfig} />
                </div>
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                    <font>Algorithms : </font>
                    &nbsp;&nbsp;
                    <Radio
                        color='primary'
                        checked={this.state.algorithm === 'singal_linkage'}
                        onChange={this.handleChange.bind(this)}
                        value="singal_linkage"
                    />
                    <font>Singal linkage</font>
                    <Radio
                        color='primary'
                        checked={this.state.algorithm === 'average_linkage'}
                        onChange={this.handleChange.bind(this)}
                        value="average_linkage"
                    />
                    <font>Average linkage</font>
                </div>
                <br />
                <div style = {{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                    <Button variant="contained" className ={classes.cssRoot}
                    onClick={this.handlePost.bind(this)}>
                        Upload
                        &nbsp;&nbsp;
                        <CloudUploadIcon />
                    </Button>
                    &nbsp;&nbsp;&nbsp;&nbsp;
                    <Link to={this.state.to} style={{ textDecoration:'none' }}>
                        <Button variant="contained" color="primary" onClick={this.submit.bind(this)}>
                            Draw &nbsp; dendrogram
                            &nbsp;&nbsp;
                            <Icon>send</Icon>
                        </Button>
                    </Link>
                </div>
                <br />
                <br />
                <br />
                <br />
            </div>
        );
    }
}

export default withStyles(styles)(Upload_profile);