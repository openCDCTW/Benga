import React from 'react';
import ReactDOM from 'react-dom';
import DropzoneComponent from 'react-dropzone-component';
import { Link } from 'react-router-dom';
import Button from '@material-ui/core/Button';
import { withStyles } from '@material-ui/core/styles';
import CloudUploadIcon from '@material-ui/icons/CloudUpload';
import Icon from '@material-ui/core/Icon';
import DeleteIcon from '@material-ui/icons/Delete';
import blue from '@material-ui/core/colors/blue';
//
import Radio from '@material-ui/core/Radio';
import SearchBar from './SearchBar.jsx';
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
            algorithm: 'single'
        };

        this.djsConfig = {
            dictDefaultMessage:"Drag cgMLST profiles here (up to 500 files)",
            dictRemoveFile:"Remove",
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

        var scheme = {};
        scheme.prof_num = fileCheck;
        scheme.linkage = this.state.algorithm;
        fetch('api/dendrogram/upload/' + window.clusteringID + '/', { 
            method:'PATCH',
            headers: new Headers({'content-type': 'application/json'}),
            body: JSON.stringify(scheme)
        });

        this.props.history.push("/dendrogram_view");
    }

    remove(){

        this.dropzone.removeAllFiles();

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

    query(){

        fetch('api/dendrogram/dendrogram/' + this.state.querybyID + '/', { method:'GET'})
        .then(response => response.json())
        .then(result => this.setState(state => ({
                png_file: result.png_file, 
                pdf_file: result.pdf_file,
                svg_file: result.svg_file, 
                emf_file: result.emf_file, 
                newick_file: result.newick_file })))
        .catch(error => alert("Data not found"));
    }

    render() {
        const config = this.componentConfig;
        const djsConfig = this.djsConfig;
        const eventHandlers = {
            init: dz => this.dropzone = dz,
        }
        const { classes } = this.props;

        if(this.state.png_file == undefined){
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
                    <div style = {{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                        <Button variant="contained" className ={classes.cssRoot}
                        onClick={this.handlePost.bind(this)}>
                            Submit
                            &nbsp;&nbsp;
                            <CloudUploadIcon />
                        </Button>
                    </div>
                    <br />
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <SearchBar
                            onChange = {(value) => this.setState({ querybyID: value })}
                            onRequestSearch={this.query.bind(this)}
                            placeholder = {"ID to get result"}
                            style = {{
                                width: '90%',
                                margin: '0 auto',
                            }}
                        />
                    </div>
                    <br />
                    <br />
                    <br />
                    <br />
                    <br />
                </div>
            );
        }else{
            return (
                <div>
                    <br />
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <img src={this.state.svg_file} />
                    </div>
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <font>Download</font>
                    </div>
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <a download href={this.state.png_file} style={{ textDecoration:'none' }}>
                            <Button variant="contained" color="default">Png</Button>
                        </a>
                        &nbsp;&nbsp;&nbsp;&nbsp;
                        <a download href={this.state.pdf_file} style={{ textDecoration:'none' }}>
                            <Button variant="contained" color="default">Pdf</Button>
                        </a>
                        &nbsp;&nbsp;&nbsp;&nbsp;
                        <a download href={this.state.svg_file} style={{ textDecoration:'none' }}>
                            <Button variant="contained" color="default">Svg</Button>
                        </a>
                        &nbsp;&nbsp;&nbsp;&nbsp;
                        <a download href={this.state.emf_file} style={{ textDecoration:'none' }}>
                            <Button variant="contained" color="default">emf</Button>
                        </a>
                        &nbsp;&nbsp;&nbsp;&nbsp;
                        <a download href={this.state.newick_file} style={{ textDecoration:'none' }}>
                            <Button variant="contained" color="default">newick</Button>
                        </a>
                    </div>
                    <br />
                    <br />
                </div>
            );
        }
    }
}

export default withStyles(styles)(Upload_profile);