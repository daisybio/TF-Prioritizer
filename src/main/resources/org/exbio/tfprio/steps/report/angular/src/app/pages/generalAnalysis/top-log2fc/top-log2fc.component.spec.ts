import { ComponentFixture, TestBed } from '@angular/core/testing';

import { TopLog2fcComponent } from './top-log2fc.component';

describe('TopLog2fcComponent', () => {
  let component: TopLog2fcComponent;
  let fixture: ComponentFixture<TopLog2fcComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ TopLog2fcComponent ]
    })
    .compileComponents();

    fixture = TestBed.createComponent(TopLog2fcComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
